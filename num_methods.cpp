#include "num_methods.h"
#include "linalg.h"
#include "graph.h"

Vec pcg(
    const function<Vec(const Vec&)>& Lx,
    Vec b,
    const function<Vec(const Vec&)>& precond,
    double eps,
    bool verbose) {
  // Solve Lx = b
  // L = BT R^{-1} B
  int max_iter = 1000000;

  int n = SZ(b);

  Vec x(n);
  Vec r = b;

  Vec z = precond(r);
  Vec p = z;
  int iter = 0;

  Vec Ap(n);
  Vec Ar(n);

  double rnorm_pre = dot(r, r);
  double rnorm = 0;
  double init_rnorm = rnorm_pre;
  double rTz_pre = dot(r, z);
  double init_rTz = dot(r, z);

  for (iter = 0; iter < max_iter; ++iter) {
    Ap = Lx(p);

    double pTAp = dot(p, Ap);
    if (verbose)
      if (iter < 10 || iter % 100 == 0) {
        Ar = Lx(r);
        double rTAr = dot(r, Ar);
        printf("iter %d |r| %g pTAp %g rTz %g rTAr %g\n",
            iter, rnorm_pre, pTAp, rTz_pre, rTAr);
      }
    //if (fabs(pTAp) <= 1e-12) {
      //if (verbose)
        //printf("PCG STOPPED due to small pTAp %g\n", pTAp);
      //break;
    //}

    double alpha = rTz_pre / pTAp;
    axpy(p, alpha, x);
    axpy(Ap, -alpha, r);
    rnorm = dot(r, r);
    if (rnorm <= eps * init_rnorm) {
      if (verbose)
        printf("CG STOPPED in %d iterations: rTr %g\n", iter, rnorm);
      break;
    }

    z = precond(r);

    double rTz = dot(r, z);

    double beta = rTz / rTz_pre;
    scal(p, beta);
    axpy(z, 1, p);

    rnorm_pre = rnorm;
    rTz_pre = rTz;
  }
  return x;
}

Vec richardson(
    const function<Vec(const Vec&)>& Lx,
    Vec b,
    const function<Vec(const Vec&)>& precond,
    double eps,
    bool verbose) {
  // Solve Lx = b
  // x' = b - L x
  // x' = P (b - L x)
  // x' = P r
  int max_iter = 1e6;
  int n = SZ(b);

  Vec r = b;
  Vec Pb = precond(b);
  Vec x(n);
  Vec Ax(n);
  Vec PAx(n);
  Vec Pr(n);
  Vec Ar(n);
  int iter = 0;

  double rnorm_pre = dot(r, r);
  double rnorm = 0;
  double init_rnorm = rnorm_pre;

  for (iter = 0; iter < max_iter; ++iter) {
    Ax = Lx(x);
    //LapVecMul(g, x, Ax);
    for (int i = 0; i < n; ++i) {
      r[i] = b[i] - Ax[i];
    }
    Pr = precond(r);

    double delta = dot(Pr, r);
    Vec APr = Lx(Pr);
    delta /= dot(Pr, APr);
    //double delta = 1.0;

    axpy(Pr, delta, x);
    // x + t Pr

    rnorm = dot(r, r);
    if (verbose)
      if (iter % 100 == 0 || iter < 10)
        printf("iter %d |r| %g |x| %g delta %g\n", iter, rnorm, dot(x, x), delta);
    if (rnorm != rnorm) {
      if (verbose)
        printf("NAN\n");
      break;
    }
    if (rnorm < eps * init_rnorm) {
      if (verbose)
        printf("CG STOPPED in %d iterations: rTr %g\n", iter, rnorm);
      break;
    }
    rnorm_pre = rnorm;
  }
  return x;
}

function<Vec(const Vec&)> buildPrecondApproxChol(const Graph& g, bool verbose) {
  int n = g.n;

  auto start = chrono::steady_clock::now();

  CholDecomp chol = approx_chol(g);

  auto end = chrono::steady_clock::now();
  if (verbose)
    printf("time spent on Approx Cholesky: %g sec\n",
        chrono::duration_cast<chrono::nanoseconds>(end - start).count() * 1e-9);

  const MatF& A = chol.A;
  MatF AT(n);
  const vector<int>& ord = chol.ord;
  for (int i = 0; i < n; ++i) {
    for (auto [j, c]: A[i])
      AT[j].push_back({i, c});
  }

  int nnz_A = 0;
  for (int i = 0; i < n; ++i)
    nnz_A += SZ(chol.A[i]);
  if (verbose)
    printf("nnz of A: %d\n", nnz_A);

  auto approx_chol_precond = [A, AT, ord](const Vec &y) {
    Vec z = chol_solve(A, AT, ord, y);
    return z;
  };

  return approx_chol_precond;
}

function<Vec(const Vec&)> buildPrecondMST(const Graph& g, bool verbose) {
  int n = g.n;

  auto weight = [&](const Edge& e) {
    return 1.0 / e.r;
  };

  auto start = chrono::steady_clock::now();

  Graph MST = mst(g, weight);

  auto end = chrono::steady_clock::now();

  if (verbose)
    printf("time spent on MST: %g sec\n",
        chrono::duration_cast<chrono::nanoseconds>(end - start).count() * 1e-9);

  if (verbose) {
    TreeQuality tq = treePrecondQuality(g, MST, weight);
    printf("tot str = %g, mx str = %g, avg str = %g, mx dila = %g\n",
        tq.tot_str, tq.mx_str, tq.tot_str / g.m, tq.mx_dila);
  }

  auto tree_precond = [MST](const Vec &y) {
    Vec z = treeSolve(MST, y);
    return z;
  };

  return tree_precond;
}

Vec lapSolve(const Graph&g, Vec b, double eps, bool verbose) {
  // Solve Lx = b
  // L = BT R^{-1} B

  //print(toLap(g));
  int n = g.n;

  double sum_eig = 0;
  double mn_r = 1e20, mx_r = 0;
  for (int i = 0; i < n; ++i)
    for (auto e: g[i]) {
      mx_r = max(mx_r, e.r);
      mn_r = min(mn_r, e.r);
      sum_eig += 1.0 / e.r;
    }

  if (verbose)
    printf("sum eigenvalues: %g min/max resis %g %g\n", sum_eig, mn_r, mx_r);

  double bT1 = 0;
  for (int i = 0; i < n; ++i)
    bT1 += b[i];

  if (verbose)
    printf("bT1 = %g\n", bT1);

  auto precond = buildPrecondApproxChol(g, verbose);
  //auto precond = buildPrecondMST(g);

  auto lapMul = [&](const Vec &x) {
    Vec Lx(n);
    LapVecMul(g, x, Lx);
    return Lx;
  };

  Vec pb = precond(b);
  //return pb;

  auto start = chrono::steady_clock::now();

  Vec x = pcg(lapMul, b, precond, eps, verbose);

  auto end = chrono::steady_clock::now();

  printf("PCG takes %g sec\n",
      chrono::duration_cast<chrono::nanoseconds>(end - start).count() * 1e-9);

  if (verbose) {
    printf("b^T P b / b^T L_inv b = %g / %g = %g\n",
        dot(b, pb), dot(b, x), dot(b, pb) / dot(b, x));
    Vec Lx = lapMul(x);
    axpy(b, -1, Lx);
    double res = dot(Lx, Lx);
    printf("|Lx - b|^2 %g\n", res);
  }
  return x;
}

Vec getElecFlow(Graph& g, const Vec& b, double eps, bool verbose) {
  int n = g.n;
  Vec x = lapSolve(g, b, eps, verbose);
  for (int u = 0; u < n; ++u)
    for (Edge& e: g[u]) {
      e.f = (x[e.u] - x[e.v]) / e.r;
    }
  return x;
}

void fixFlow(Graph& g, const Vec b) {
  int n = g.n;
  Vec remain = b;
  for (int u = 0; u < n; ++u)
    for (auto e: g[u]) {
      remain[u] -= e.f;
    }

  vector<bool> vis(n);
  vector<int> mom(n);
  vector<double> demand_subtree(n);

  function<void(int)> bfs = [&](int s) {
    //printf("dfs %d %d %g\n", u, p, demand);
    vector<int> qu;
    int ql = 0;
    vis[s] = 1;
    mom[s] = s;
    qu.push_back(s);
    while (ql < SZ(qu)) {
      int u = qu[ql++];
      for (Edge& e: g[u]) if (!vis[e.v]) {
        vis[e.v] = 1;
        mom[e.v] = u;
        qu.push_back(e.v);
      }
    }
    reverse(qu.begin(), qu.end());
    for (int u: qu) {
      demand_subtree[u] += remain[u];
      for (Edge &e: g[u]) if (e.v == mom[u]) {
        e.f            += demand_subtree[u];
        g[e.v][e.ri].f -= demand_subtree[u];
        demand_subtree[e.v] += demand_subtree[u];
      }
    }
  };

  for (int u = 0; u < n; ++u) {
    if (vis[u])
      continue;
    bfs(u);
  }
}

/*
 * Compute min_{d: B^T d = 0, gT d = c} d^T R d
 * What if no solution??
 */
Vec getGradientCirc(Graph& g, double c, double eps) {
  int n = g.n;

  // uu = BT R^{-1} g
  Vec uu(n);
  for (int u = 0; u < n; ++u)
    for (auto e: g[u]) {
      uu[u] += e.g / e.r;
    }

  //auto Linv = buildPrecondApproxChol(g, true);
  //Vec Linv_uu = Linv(uu);
  Vec Linv_uu = lapSolve(g, uu, 1e-24, false);

  // vv = (c / (u^T L^{-1} u - a)) L^{-1} u
  Vec vv = Linv_uu;
  Vec Lvv(g.n);
  LapVecMul(g, vv, Lvv);
  axpy(uu, -1, Lvv);
  printf("|L vv - uu|_2^2 = %g\n", dot(Lvv, Lvv));

  // f = R^{-1} (B vv + y g)
  for (int u = 0; u < n; ++u) {
    for (Edge &e: g[u]) {
      e.f = (vv[e.u] - vv[e.v] - e.g) / e.r;
    }
  }

  // Fix the flow via tree

  fixFlow(g, Vec(n, 0));

  double gTf = 0;
  for (int u = 0; u < n; ++u) {
    for (auto e: g[u]) if (e.u < e.v) {
      gTf += e.g * e.f;
    }
  }
  for (int u = 0; u < n; ++u) {
    for (Edge& e: g[u]) {
      e.f *= c / gTf;
    }
  }

  return vv;
}

/*
 * Compute sum_e r_e f_e^2
 */
double l2Cost(const Graph& g) {
  int n = g.n;
  double res = 0;
  for (int u = 0; u < n; ++u)
    for (const Edge& e: g[u]) if (e.u < e.v) {
      res += e.r * (e.f * e.f);
    }
  return res;
}

/*
 * Compute sum_e s_e^p f_e^p
 */
double pnormCost(const Graph& g, int p) {
  int n = g.n;
  double res = 0;
  for (int u = 0; u < n; ++u)
    for (const Edge& e: g[u]) if (e.u < e.v) {
      res += pow(e.s * fabs(e.f), p);
    }
  return res;
}


/*
 * Compute min_{d: B^T d = 0, gT S d = i/2} d^T S^T (R + s I) S d
 * where
 * R = |S f|^{p - 2}
 * g = p R S f
 * s = 0.5 i^{(p - 2) / p} m^{-(p - 2) / p}
 */
Graph findDelta(Graph g, double i, int p, double eps=1e-9) {
  int n = g.n;
  int m = 0;
  for (int u = 0; u < n; ++u)
    m += SZ(g[u]);
  m /= 2;

  Graph g2 = g;

  double s = 0.5 * pow(i, 1.0 * (p - 2) / p) / pow(m, 1.0 * (p - 2) / p);

  for (int u = 0; u < n; ++u) {
    for (Edge &e: g2[u]) {
      double r = pow(fabs(e.s * e.f), p - 2);
      e.r = e.s * (r + s) * e.s;
      e.g = e.s * p * r * e.s * e.f;
    }
  }

  getGradientCirc(g2, i / 2.0, eps);

  return g2;
}

double lineSearch(Graph g, Graph delta, int p) {
  int n = g.n;
  auto grad = [&](double t) {
    /*
     * (d / dt) sum_e s_e^p (f_e - t delta_e)^p
     * = sum_e s_e^p p (f_e - t delta_e)^{p - 1} (-delta_e)
     */
    double res = 0;

    for (int u = 0; u < n; ++u) {
      assert(SZ(g[u]) == SZ(delta[u]));
      for (int i = 0; i < SZ(g[u]); ++i) {
        Edge e = g[u][i];
        if (e.u > e.v)
          continue;
        double delta_e = delta[u][i].f;
        double nf_e = e.f - t * delta_e;
        double scale = pow(e.s, p);
        double grad_e = scale * p * pow(fabs(nf_e), p - 2) * nf_e * (-delta_e);
        res += grad_e;
      }
    }

    return res;
  };

  double L = -1, R = 1;
  while (grad(R) < 0) {
    L = R;
    R = 2 * R;
  }
  while (grad(L) > 0) {
    R = L;
    L = 2 * L;
  }

  for (int iter = 0; iter < 50; ++iter) {
    double M = (L + R) / 2;
    if (grad(M) < 0)
      L = M;
    else
      R = M;
  }
  return L;
}

void takeStep(Graph& g, const Graph& delta, double alpha) {
  int n = g.n;
  for (int u = 0; u < n; ++u) {
    assert(SZ(g[u]) == SZ(delta[u]));
    for (int i = 0; i < SZ(g[u]); ++i) {
      g[u][i].f -= alpha * delta[u][i].f;
    }
  }
}

double evalResidual(const Graph &g, const Graph &delta, const double alpha, const int p) {
  int n = g.n;

  double res = 0;
  for (int u = 0; u < n; ++u) {
    assert(SZ(g[u]) == SZ(delta[u]));
    for (int i = 0; i < SZ(g[u]); ++i) {
      Edge e = g[u][i];
      double df = delta[u][i].f * alpha;
      if (e.u > e.v)
        continue;

      double r = pow(e.s * fabs(e.f), p - 2);
      double eg = p * r * e.s * e.f;

      double tmp = eg * e.s * df;
      tmp -= 2.0 * p * p * (df * e.s * r * e.s * df);
      tmp -= pow(p, p) * pow(e.s * df, p);
    }
  }
  return res;
}

/*
 * delta_e.s = g_e.s
 * delta_e.g = g_e.s * p * (g_e.s * g_e.f)^{p - 1}
 * delta_e.r = g_e.s^2 * (s + (g_e.s * g_e.f)^{p - 2})
 */
bool progressCheck(const Graph &g, const Graph &delta, const double i, const int p) {
  int n = g.n;
  int m = 0;
  for (int u = 0; u < n; ++u)
    m += SZ(g[u]);
  m /= 2;

  double lambda = 16 * p;

  double delta_pnorm_val = pnormCost(delta, p);
  double delta_l2_val = l2Cost(delta);
  double k = pow(p, p) * delta_pnorm_val / (2.0 * p * p * delta_l2_val);

  double alpha_0 = min(1.0 / (16 * lambda), 1.0 / pow(16.0 * lambda * k, 1.0 / (p - 1)));

  if (delta_l2_val > lambda * i / (p * p))
    return 0;

  if (evalResidual(g, delta, alpha_0, p) < alpha_0 * i / 4.0)
    return 0;

  return 1;
}

double grad(const Graph& g, double i, int p) {
  int n = g.n;
  int m = 0;
  for (int u = 0; u < n; ++u)
    m += SZ(g[u]);
  m /= 2;

  Graph g2 = g;
  double s = 0.5 * pow(i, 1.0 * (p - 2) / p) / pow(m, 1.0 * (p - 2) / p);

  double gTWg = 0;
  for (int u = 0; u < n; ++u) {
    for (Edge &e: g2[u]) {
      double r = pow(fabs(e.s * e.f), p - 2);
      e.r = e.s * (r + s) * e.s;
      e.g = e.s * p * r * e.s * e.f;
      if (e.u < e.v)
        gTWg += e.g * e.g / e.r;
    }
  }

  return gTWg;
}

double violation(Graph g, Vec b) {
  int n = g.n;
  double res = 0;
  for (int u = 0; u < n; ++u) {
    double uf = 0;
    for (auto e: g[u])
      uf += e.f;
    res = max(res, fabs(b[u] - uf));
  }
  return res;
}

/*
 * Compute min_f sum_e s_e^p |f_e|^p = |s f|_p^p s.t. B^Tf = b
 */
Graph pIRLS(Graph g, Vec b, int p, double eps) {
  int n = g.n;

  int m = 0;
  for (int u = 0; u < n; ++u)
    m += SZ(g[u]);
  m /= 2;

  Graph g0 = g;
  for (int u = 0; u < n; ++u)
    for (Edge& e: g0[u]) {
      e.r = e.s * e.s;
    }
  getElecFlow(g0, b, 1e-22, false);
  fixFlow(g0, b);

  double cur_pnorm = pnormCost(g0, p);

  double i = cur_pnorm / (16. * p);

  printf("|S f2|_p^p = %g, i = %g\n", cur_pnorm, i);

  double kappa = 2.0 * eps / (16. * p * (1 + eps));


  Graph cur_g = g0;
  for (int u = 0; u < n; ++u)
    for (Edge& e: cur_g[u]) {
      e.r = 0;
    }

  double init_grad = grad(cur_g, i / 2, p);

  if (init_grad < pow(eps, p))
    return cur_g;

  for (int iter = 0; cur_pnorm * kappa < i; ++iter) {
    double cur_grad = grad(cur_g, i / 2, p);
    double max_vio = violation(cur_g, b);

    printf("iter %d |S f|_p^p = %g |S f|_p = %g i = %g, grad = %g, vio = %g\n", iter, cur_pnorm, pow(cur_pnorm, 1.0 / p), i, cur_grad, max_vio);

    if (cur_grad <= pow(eps, p) * init_grad) {
      printf("grad too small %g\n", cur_grad);
      break;
    }
    Graph delta = findDelta(cur_g, i / 2, p);

    double alpha = lineSearch(cur_g, delta, p);

    takeStep(cur_g, delta, alpha);

    cur_pnorm = pnormCost(cur_g, p);

    if (!progressCheck(cur_g, delta, i, p))
      i /= 2;
  }

  return cur_g;
}
