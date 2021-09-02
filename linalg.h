#ifndef LINALG
#define LINALG

#include "common.h"
#include "graph_instance.h"

void LapVecMul(const Graph &g, const Vec& x, Vec& y);

double dot(const Vec &a, const Vec &b);

void axpy(const Vec &x, double a, Vec& y);

void scal(Vec &a, double t);

Mat toLap(const Graph &g);

Mat recover(Mat A);

void print(const Mat& A);

void printF(const MatF& A);

Vec matMul(const MatF& A, const Vec& x);

struct CholDecomp {
  //Vec D;
  MatF A;
  vector<int> ord;
};

Vec chol_solve(const MatF& A, const MatF& AT, const vector<int>& ord, const Vec& b);

/*
 * Output ord, A s.t.
 * A^T A \approx L
 */
CholDecomp approx_chol(Graph);

#endif
