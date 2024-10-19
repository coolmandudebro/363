#ifndef LAB_5_HEADER_H
#define LAB_5_HEADER_H
#include <iostream>
#include "C:/Users/gerge/CLionProjects/eigen-3.4.0/Eigen/Eigen"
#include <cmath>
#include <fstream>
#include <stdio.h>
using namespace std;
using namespace Eigen;
using Eigen::MatrixXd;
//Setup

// Observations Matrix
void read(MatrixXd& A,  char* file, int rows, int cols);

// Scientific Wild Ass Guess
void swag(MatrixXd& x, MatrixXd l);

// Better Printf
void print(MatrixXd M, char* name, int units);

// A Design Matrix
MatrixXd designA(int m, int n, MatrixXd l, MatrixXd x);

// B Design Matrix
MatrixXd designB(int m, int n, MatrixXd l, MatrixXd x);

// Functional Model
MatrixXd funcModel(MatrixXd l, MatrixXd x);

// Adjustment Loop
void implicit(MatrixXd& delta, MatrixXd& v, double thresh, MatrixXd& A, MatrixXd& B, MatrixXd& x, MatrixXd l,
              MatrixXd& ln, MatrixXd& w, MatrixXd& N, MatrixXd& M, MatrixXd P, int n, int u);

// Correlation Matrix
MatrixXd corr(MatrixXd C);

// Latex Matrix Print
void square_matrix(MatrixXd M, char* name1, char* name2, bool hat, int units, FILE *output);

// Latex Table Print
void vectorf(MatrixXd M, char* name1, char* name2, bool hat, int units, FILE *output);

#endif //LAB_5_HEADER_H
