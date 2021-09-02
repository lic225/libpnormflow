#ifndef GRAPH_INSTANCE
#define GRAPH_INSTANCE
#include "common.h"

struct Edge {
  int u, v;
  int ri;
  double g, r, s; // gradient, resistance, and p norm weight
  double f;
};

typedef function<double(const Edge&)> EFunc;

/*
 * Undirected weighted graph
 */
struct Graph {
  int n, m;
  vector<vector<Edge>> G;
  void init(int _n, int _m) {
    n = _n;
    m = _m;
    G = vector<vector<Edge>>(n);
  }
  void addEdge(int u, int v, double g=0, double r=0, double s=1) {
    G[u].push_back({u, v, (int)G[v].size(), g, r, s, 0});
    G[v].push_back({v, u, (int)G[u].size() - 1, g, r, s, 0});
  }
  vector<Edge>& operator[](int u) { return G[u]; }
  const vector<Edge>& operator[](int u) const { return G[u]; }
};
#endif
