#include "common.h"
#include "graph_instance.h"
#include "linalg.h"
#include "graph.h"
#include "num_methods.h"
//#include "maxflow_comb.h"

//Vec pIRLS(const Graph &g) {
//}

int main(int argc, char** argv) {
  Graph g;
  int n, m;
  Vec b;
  bool RTE = argc > 1 && strcmp(argv[1], "RTE") == 0;

  int p;

  if (RTE) {
    scanf("%d", &n);
    g.init(n, 0);
    m = 0;

    for (int i = 0; i < n; ++i) {
      int deg;
      scanf("%d", &deg);
      m += deg;
      int cur = i + 1;
      for (int j = 0; j < deg; ++j) {
        int v;
        double w;
        scanf("%d%lf", &v, &w);
        cur += v;
        g.addEdge(i, cur-1, 0, 1.0 / w, 1.0 / sqrt(w));
      }
    }
    b = Vec(n);
    for (int i = 0; i < n; ++i) {
      scanf("%lf", &b[i]);
    }
    g.m = m;

  } else {
    scanf("%d%d%d", &n, &n, &m);

    g.init(n, m);

    for (int i = 0; i < m; ++i) {
      int u, v;
      scanf("%d%d", &u, &v);
      //++u, ++v;
      g.addEdge(u, v, 0, 1, 1);
    }

    b = Vec(n);
    b[0] = 1;
    b[n - 1] = -1;
  }

  printf("# of edges: %d\n", m);
  auto [n_comp, comp_ids] = components(g);
  printf("# of components : %d\n", n_comp);


  //Vec x = lapSolveTree(g, b, 0.1);
  Vec x = getElecFlow(g, b, 1e-22, true);
  fixFlow(g, b);
  //Vec x = getElecFlow(g, b, 1e-22, true);
  double ans = dot(x, b);
  printf("%.12f\n", ans);

  double l2_cong = 0;
  double l2_energy = 0;
  for (int u = 0; u < n; ++u)
    for (auto e: g[u]) if (e.u < e.v) {
      l2_cong = max(l2_cong, fabs(e.f) * sqrt(e.r));
      l2_energy += e.f * e.f * e.r;
    }
  printf("L2 cong = %g\n", l2_cong);
  //printf("L2 flow value = %g\n", 1.0 / l2_cong);
  printf("L2 energy = %g\n", l2_energy);

  if (m < 1e5) {
    for (int u = 0; u < n; ++u)
      for (Edge& e: g[u])
        e.f = 0;
    Graph g_inf;
    int src = n;
    int snk = n + 1;
    double L = 0, R = 1e9;
    for (int iter = 0; iter < 100; ++iter) {
      double M = (L + R) / 2;
      g_inf = g;
      g_inf.n += 2;
      g_inf.G.push_back({});
      g_inf.G.push_back({});
      double total_in = 0;
      for (int i = 0; i < n; ++i)
        if (b[i] > 0) {
          total_in += M * b[i];
          double cap = M * b[i];
          g_inf.addEdge(src, i, 0, 1.0 / (cap * cap), 0);
        } else if (b[i] < 0) {
          double cap = -M * b[i];
          g_inf.addEdge(i, snk, 0, 1.0 / (cap * cap), 0);
        }
      double mxf = maxflow_isap(g_inf, src, snk, [&](const Edge &e) {return 1.0 / sqrt(e.r);});
      //printf("try M=%g with total in %g, mxf=%.9f\n", M, total_in, mxf);
      if (mxf < total_in - 1e-9) {
        R = M;
      } else {
        L = M;
      }
    }

    double linf_cong = 0;
    double linf_energy = 0;
    for (int u = 0; u < n; ++u)
      for (auto e: g_inf[u]) if (e.v != src && e.v != snk) {
        double f = e.f / L;
        linf_cong = max(linf_cong, fabs(f) * sqrt(e.r));
        linf_energy += f * f * e.r;
      }
    printf("M = %g\n", L);
    printf("Linf cong = %g\n", linf_cong);
    printf("Linf energy = %g\n", linf_energy);
  }

  p = max(2, (int)ceil(sqrt(log2(m))));
  //p = ceil(log2(m) / 0.5);
  printf("p=%d\n", p);
  for (int u = 0; u < n; ++u)
    for (Edge &e: g[u])
      e.g = e.r = 0;
  Graph g_p = pIRLS(g, b, p, 1e-6);
  double lp_energy = pnormCost(g_p, p);
  printf("|S fp|_p^p = %g\n", lp_energy);
  printf("|S fp|_p = %g\n", pow(lp_energy, 1.0 / p));
  
  double lp_cong = 0;
  for (int u = 0; u < n; ++u)
    for (auto e: g_p[u])
      lp_cong = max(lp_cong, fabs(e.f * e.s));
  printf("|S fp|_inf = %g\n", lp_cong);
}
