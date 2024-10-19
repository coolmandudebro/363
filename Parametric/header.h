#ifndef LAB_4_HEADER_H
#define LAB_4_HEADER_H

#include <iostream>
#include "C:/Users/gerge/CLionProjects/eigen-3.4.0/Eigen/Eigen"
#include <cmath>
#include <fstream>
#include <stdio.h>
using namespace std;
using namespace Eigen;
using Eigen::MatrixXd;
//Setup

void Lread(MatrixXd& A,  char* file);

void Uread(MatrixXd& A,  char* file);

void print(MatrixXd M, char* name, int units);

void parametric(MatrixXd& delta, double thresh, MatrixXd& A, MatrixXd& x, MatrixXd c,
                MatrixXd l, MatrixXd& w, MatrixXd& N, MatrixXd P, int n, int u);

MatrixXd corr(MatrixXd C);

void square_matrix(MatrixXd M, char* name1, char* name2, bool hat, int units, char* file);

void vector(MatrixXd M, char* name1, char* name2, int units, char* file);





#endif //LAB_4_HEADER_H
