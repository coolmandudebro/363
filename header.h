#ifndef LAB_3_HEADER_H
#define LAB_3_HEADER_H
#include <iostream>
#include "C:/Users/gerge/CLionProjects/eigen-3.4.0/Eigen/Eigen"
#include <cmath>
#include <fstream>
#include <stdio.h>
using namespace std;
using namespace Eigen;
using Eigen::MatrixXd;
//Setup

void Bread(MatrixXd& A,  char* file);

void Cread(MatrixXd& A,  char* file);

void Lread(MatrixXd& A,  char* file);

void Jread(MatrixXd& A,  char* file);

void print(MatrixXd M, char* name, int units);

void experimentalprint(MatrixXd M, char* name, int units);

MatrixXd corr(MatrixXd C);

void fprint(MatrixXd M, int units, char* file);

#endif //LAB_3_HEADER_H
