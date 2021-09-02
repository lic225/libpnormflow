#ifndef GRAPH
#define GRAPH

#include "common.h"
#include "graph_instance.h"


vector<double> shortestPath(const Graph &g, int source, const EFunc& length);

pair<int, vector<int>> components(const Graph &g);

double maxflow_isap(Graph &g, int source, int sink, const EFunc& cap);

Graph mst(const Graph &g, const EFunc& length);

struct TreeQuality {
  double tot_str;
  double mx_str;
  double mx_dila;
};

TreeQuality treePrecondQuality(
    const Graph &g, const Graph &T, const EFunc& length);

Vec treeSolve(const Graph &T, const Vec &y);

#endif
