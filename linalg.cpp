#include "linalg.h"

void LapVecMul(const Graph &g, const Vec& x, Vec& y) {
  // y = Lx
  y.resize(g.n);
  for (int i = 0; i < g.n; ++i) {
    y[i] = 0;
    for (const Edge &e: g[i]) {
      y[i] += (x[i] - x[e.v]) / e.r;
    }
  }
}

double dot(const Vec &a, const Vec &b) {
  assert(a.size() == b.size());
  double res = 0;
  for (size_t i = 0; i < a.size(); ++i)
    res += a[i] * b[i];
  return res;
}

void axpy(const Vec &x, double a, Vec& y) {
  // y = a * x + y
  assert(x.size() == y.size());
  for (size_t i = 0; i < y.size(); ++i)
    y[i] += a * x[i];
}

void scal(Vec &a, double t) {
  // a =  t * a
  for (size_t i = 0; i < a.size(); ++i)
    a[i] = t * a[i];
}

Mat toLap(const Graph &g) {
  Mat m(g.n);
  for (int i = 0; i < g.n; ++i)
    for (const Edge &e: g[i]) {
      m[i][e.v] -= 1.0 / e.r;
      m[i][i] += 1.0 / e.r;
    }
  return m;
}

Mat recover(Mat A) {
  // A^T A = sum_u A[u]^T A[u]
  Mat res(SZ(A));
  for (int i = 0; i < SZ(A); ++i) {
    for (auto [j, c1]: A[i])
      for (auto [k, c2]: A[i]) {
        res[j][k] += c1 * c2;
      }
  }
  return res;
}

void print(const Mat& A) {
  for (int i = 0; i < SZ(A); ++i) {
    printf("%d: ", i);
    for (auto [v, c]: A[i])
      printf("%d:%g ", v, c);
    puts("");
  }
}

void printF(const MatF& A) {
  for (int i = 0; i < SZ(A); ++i) {
    printf("%d: ", i);
    for (auto [v, c]: A[i])
      printf("%d:%g ", v, c);
    puts("");
  }
}

Vec matMul(const MatF& A, const Vec& x) {
  Vec res(SZ(A));
  for (int i = 0; i < SZ(A); ++i)
    for (auto [j, c]: A[i])
      res[i] += c * x[j];
  return res;
}

/*
 * Output ord, A s.t.
 * A^T A \approx L
 */
CholDecomp approx_chol(Graph g) {
  int n = g.n;
  int m = 0;
  for (int i = 0; i < n; ++i)
    m += SZ(g[i]);
  m /= 2;

  mt19937 prg(514514);
  vector<int> ord(n);
  for (int i = 0; i < n; ++i)
    ord[i] = i;

  random_shuffle(ord.begin(), ord.end());
  //sort(ord.begin(), ord.end(), [&](int u, int v){return SZ(g[u]) < SZ(g[v]);});
  vector<int> inv_ord(n);
  for (int i = 0; i < n; ++i)
    inv_ord[ord[i]] = i;

  int k = max(12, 1000000 / m);
  //int K = 5;
  //g = split(g, k);
  vector<vector<pair<int,double>>> G(n);

  for (int u = 0; u < n; ++u) {
    for (Edge e: g[u]) {
      int v = e.u == u? e.v: e.u;
      if (inv_ord[u] < inv_ord[v]) {
        double er = e.r;
        for (int j = 0; j < k; ++j) {
          G[u].emplace_back(v, k * er);
        }
      }
    }
  }


  Vec D(n);
  MatF A(n);
  //return CholDecomp{A, ord};

  vector<double> c_sum;
  map<int, double> eu;

  for (int i = 0; i < n - 1; ++i) {
    int u = ord[i];

    if (SZ(G[u]) == 0) continue;

    D[u] = 0;
    for (auto [v, r]: G[u]) {
      D[u] += 1.0 / r;
    }
    double sqrtD = sqrt(D[u]);

    eu.clear();
    for (auto [v, r]: G[u])
      eu[v] -= (1.0 / r / sqrtD);
    eu[u] = sqrtD;
    for (auto [v, r]: eu)
      A[u].emplace_back(v, r);

    //sort(G[u].begin(), G[u].end());
    //bool pushAu = 0;
    //for (int j = 0; j < SZ(G[u]);) {
      //int v = G[u][j].first;
      //double tmp = 0;
      //for (; G[u][j].first == v && j < SZ(G[u]); ++j) {
        //double r = G[u][j].second;
        //tmp += (1.0 / r) / sqrtD;
      //}
      //if (!pushAu && u < v) {
        //pushAu = 1;
        //A[u].emplace_back(u, sqrtD);
      //}
      //A[u].emplace_back(v, -tmp);
    //}

    c_sum.clear();
    double cur = 0;
    for (auto [v, r]: G[u]) {
      cur += 1. / r;
      c_sum.push_back(cur);
    }

    uniform_real_distribution<double> real_dist(0.0, c_sum.back());
    uniform_int_distribution<int> range_dist(0, SZ(G[u]) - 1);
    for (int j = 0; j < SZ(G[u]); ++j) {
      double x = real_dist(prg);
      int pos_wgt = lower_bound(c_sum.begin(), c_sum.end(), x) - c_sum.begin();
      //int pos_wgt = range_dist(prg);
      int pos_uni = range_dist(prg);
      //int pos_uni = j;
      //int pos_wgt = rand() % SZ(G[u]);
      //int pos_uni = rand() % SZ(G[u]);

      auto [v1, r1] = G[u][pos_wgt];
      auto [v2, r2] = G[u][pos_uni];

      if (v1 != v2) {
        double nr = r1 + r2;
        if (inv_ord[v1] < inv_ord[v2]) {
          G[v1].emplace_back(v2, nr);
        } else {
          G[v2].emplace_back(v1, nr);
        }
      }
    }
  }
  return CholDecomp{A, ord};
}

Vec chol_solve(const MatF& A, const MatF& AT, const vector<int>& ord, const Vec& b) {
  // Solve A^T A x = b, A: n x n upper tri
  int n = SZ(ord);


  //printF(A);
  // Solve A^T z = b
  Vec z(n);
  //printF(AT);
  for (int i = 0; i < n; ++i) {
    int u = ord[i];
    double cur = b[u];
    double ATuu = 0;
    for (auto [v, c]: AT[u]) {
      if (v != u)
        cur -= c * z[v];
      else
        ATuu = c;
    }
    if (ATuu != 0)
      z[u] = cur / ATuu;
  }
  // Solve Ax = z
  Vec x(n);
  for (int i = n-1; i >= 0; --i) {
    int u = ord[i];
    double cur = z[u];
    double Auu = 0;
    for (auto [v, c]: A[u]) {
      if (v != u)
        cur -= c * x[v];
      else
        Auu = c;
    }
    if (Auu != 0)
      x[u] = cur / Auu;
  }

  //scal(x, 0.1);

  return x;
}

