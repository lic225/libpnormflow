#include "graph.h"


vector<double> shortestPath(const Graph &g, int source, const function<double(const Edge&)>& length) {
  vector<double> dist(g.n, 1e9);
  set<pair<double, int>> pq{{0, source}};
  while (pq.size()) {
    auto [d, u] = *pq.begin();
    pq.erase(pq.begin());
    if (d < dist[u]) {
      dist[u] = d;
      for (const Edge& e: g[u]) {
        pq.insert({d + length(e), e.v});
      }
    }
  }
  return dist;
}

pair<int, vector<int>> components(const Graph &g) {
  vector<bool> vis(g.n);
  vector<int> comp_id(g.n);
  int n_comp = 0;
  function<void(int, int)> bfs = [&](int s, int id) {
    //printf("dfs %d %d\n", u, dep);
    vis[s] = 1;
    vector<int> qu{s};
    size_t ql = 0;
    while (ql < qu.size()) {
      int u = qu[ql++];
      comp_id[u] = id;
      for (auto& e: g[u]) {
        if (!vis[e.v]) {
          vis[e.v] = 1;
          qu.push_back(e.v);
        }
      }
    }
  };

  for (int i = 0; i < g.n; ++i) {
    if (!vis[i]) {
      //printf("dfs %d\n", i);
      bfs(i, ++n_comp);
    }
  }
  return {n_comp, comp_id};
}

double maxflow_isap(Graph &g, int source, int sink, const EFunc& cap) {
  vector<int> d(g.n);
  vector<size_t> iter(g.n);
  vector<int> gap(g.n);

  function<double(int, double)> push = [&](int u, double flow) {
    if (u == sink) {
      return flow;
    }
    for (size_t &i = iter[u]; i < g[u].size(); ++i) {
      Edge &e = g[u][i];
      double res = cap(e) - e.f;
      if (d[e.v] + 1 == d[u] && res > 0) {
        double f = push(e.v, min(flow, res));
        if (f > 0) {
          e.f += f;
          g[e.v][e.ri].f -= f;
          return f;
        }
      }
    }
    if( (--gap[d[u]]) == 0) d[source] = g.n;
    else {
      d[u]++;
      iter[u] = 0;
      ++gap[d[u]];
    }
    return 0.0;
  };

  double mxf = 0;
  gap[0] = g.n;
  int n_iter = 0;
  while (d[source] < g.n) {
    double cf = push(source, 1e9);
    n_iter++;
    mxf += cf;
  }
  //printf("# of iter %d\n", n_iter);
  return mxf;
}

struct Union {
  vector<int> mom;
  Union(int n) {
    init(n);
  }
  void init(int n) {
    mom = vector<int>(n);
    for (int i = 0; i < n; ++i)
      mom[i] = i;
  }
  int find(int x) {
    if (mom[x] == x)
      return x;
    else
      return mom[x] = find(mom[x]);
  }
  void merge(int a, int b) {
    mom[find(a)] = find(b);
  }
};

Graph mst(const Graph &g, const EFunc& length) {
  int n = g.n;
  Union uni(n);
  Graph T;
  vector<Edge> es;
  printf("mst n=%d\n", n);
  for (int u = 0; u < n; ++u)
    for (auto e: g[u]) if (e.u < e.v)
      es.push_back(e);

  sort(es.begin(), es.end(), [&] (const Edge &a, const Edge &b) {
    return length(a) < length(b);
  });

  double tot = 0;

  T.init(n, 0);
  for (auto e: es) {
    if (uni.find(e.u) != uni.find(e.v)) {
      //printf("merge %d %d\n", e.u, e.v);
      uni.merge(e.u, e.v);
      tot += length(e);
      T.addEdge(e.u, e.v, e.g, e.r, e.s);
    }
  }
  printf("mst cost %g\n", tot);

  return T;
}

/*
 * Compute the total and the max stretch of G w.r.t. T
 */
TreeQuality treePrecondQuality(
    const Graph &g, const Graph &T, const EFunc& length) {
  int n = g.n;
  int lg = 0;
  while ((1 << lg) < n) ++lg;
  printf("n=%d lg=%d\n", n, lg);

  vector<vector<int>> anc(n, vector<int>(lg, -1));
  vector<int> lvl(n);
  vector<double> dist(n);
  vector<int> pre_ord;

  function<void(int, int)> dfs = [&](int u, int p) {
    anc[u][0] = p;
    pre_ord.push_back(u);
    for (auto e: T[u]) if (e.v != p) {
      dist[e.v] = dist[u] + length(e);
      lvl[e.v] = lvl[u] + 1;
      dfs(e.v, u);
    }
  };

  auto lca = [&](int u, int v) -> int {
    if (lvl[u] > lvl[v]) swap(u, v);
    int jmp = lvl[v] - lvl[u];
    for (int j = 0; j < lg; ++j)
      if ((jmp >> j) & 1)
        v = anc[v][j];
    if (u == v)
      return u;

    for (int j = lg-1; j >= 0; --j)
      if (anc[u][j] != anc[v][j]) {
        u = anc[u][j];
        v = anc[v][j];
      }
    if (u == v)
      return u;
    else
      return -1;
  };

  auto d = [&](int u, int v) -> double {
    int a = lca(u, v);
    if (a == -1)
      return -1;
    double res = dist[u] + dist[v] - 2 * dist[a];
    return res;
  };

  for (int i = 0; i < n; ++i)
    if (anc[i][0] == -1) {
      dfs(i, i);
    }

  for (int j = 0; j + 1 < lg; ++j)
    for (int i = 0; i < n; ++i)
      anc[i][j+1] = anc[anc[i][j]][j];

  vector<double> str_sum(n);

  double tot = 0, mx = 0;
  for (int u = 0; u < n; ++u) {
    for (auto e: g[u]) if (e.u < e.v) {
      double str_e = d(e.u, e.v) / length(e);
      mx = max(str_e, mx);
      tot += str_e;
      int a = lca(e.u, e.v);
      str_sum[e.u] += str_e;
      str_sum[e.v] += str_e;
      str_sum[a] -= 2 * str_e;
    }
  }

  double mx_dila = 0;
  for (int i = n - 1; i >= 0; --i) {
    int u = pre_ord[i];
    str_sum[anc[u][0]] += str_sum[u];
    mx_dila = max(mx_dila, str_sum[u]);
  }

  return TreeQuality{
    .tot_str = tot,
    .mx_str = mx,
    .mx_dila = mx_dila
  };
}

Vec treeSolve(const Graph &T, const Vec& y) {
  int n = T.n;
  Vec x(n);

  vector<int> deg(n);
  for (int u = 0; u < n; ++u)
    deg[u] = SZ(T[u]);

  vector<int> qu;
  int ql = 0;
  for (int u = 0; u < n; ++u)
    if (deg[u] == 1)
      qu.push_back(u);

  vector<bool> vis(n);
  vector<int> prv(n);
  Vec ty = y;
  Vec c(n);

  while (ql < SZ(qu)) {
    int u = qu[ql++];
    vis[u] = 1;
    prv[u] = u;
    for (auto e: T[u]) if (!vis[e.v]) {
      prv[u] = e.v;
      c[u] = ty[u] * e.r;
      deg[e.v]--;
      ty[e.v] += ty[u];
      if (deg[e.v] == 1) {
        qu.push_back(e.v);
      }
      break;
    }
  }

  reverse(qu.begin(), qu.end());

  for (int i = 0; i < SZ(qu); ++i) {
    int u = qu[i];
    x[u] = x[prv[u]] + c[u];
  }

  return x;
}
