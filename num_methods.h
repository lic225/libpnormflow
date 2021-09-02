#ifndef NUM_METHODS
#define NUM_METHODS

#include "common.h"
#include "graph_instance.h"

Vec pcg(const function<Vec(const Vec&)>& Lx, Vec b, const function<Vec(const Vec&)>& precond, double eps=1e-6, bool verbose=true);

Vec richardson(const function<Vec(const Vec&)>& Lx, Vec b, const function<Vec(const Vec&)>& precond, double eps=1e-9, bool verbose=true);

// min_f sum_e r_e |f_e|^2 s.t. B^Tf = b
Vec lapSolve(const Graph&g, Vec b, double eps=1e-2, bool verbose=true);

Vec getElecFlow(Graph& g, const Vec& b, double eps=1e-9, bool verbose=true);
void fixFlow(Graph& g, const Vec b);

// sum_e s_e^p f_e^p
double pnormCost(const Graph& g, int p = 10);

// min_f sum_e s_e^p |f_e|^p s.t. B^Tf = b
Graph pIRLS(Graph g, Vec b, int p = 10, double eps=1e-2);

//Vec lapSolveTree(const Graph&g, Vec b, double eps=1e-2);

#endif
