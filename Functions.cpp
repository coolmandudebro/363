#include "header.h"

// Design Matrix
void Bread(MatrixXd& A,  char* file){

    //Populate eigen matrix
    A.resize(9,20);

    //Pointer magic
    char fileSpec[strlen(file)+1];
    snprintf( fileSpec, sizeof( fileSpec ), "%s", file );
    FILE *out = fopen( fileSpec, "r" );


    //Read and store formatted data to Matrics
    for (int i = 0; i < 9; ++i) {
        fscanf(out,"%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,\n",
               &A(i,0), &A(i,1), &A(i,2), &A(i,3), &A(i,4),
               &A(i,5), &A(i,6), &A(i,7), &A(i,8), &A(i,9),
               &A(i,10), &A(i,11), &A(i,12), &A(i,13), &A(i,14),
               &A(i,15), &A(i,16), &A(i,17), &A(i,18), &A(i,19));
    }

    fclose(out);
}

// Variance of Observations Matrix
void Cread(MatrixXd& A,  char* file){

    //Populate eigen matrix
    A.setZero(20,20);

    //Pointer magic
    char fileSpec[strlen(file)+1];
    snprintf( fileSpec, sizeof( fileSpec ), "%s", file );
    FILE *out = fopen( fileSpec, "r" );


    //Read and store formatted data to Matrics
    for (int i = 0; i < 20; ++i) {
        fscanf(out,"%lf\n",&A(i,i));
        // Convert cm to m^2
        A(i,i) = A(i,i)*A(i,i)/10000;
    }

    fclose(out);
}

// Observations Matrix
void Lread(MatrixXd& A,  char* file){

    //Populate eigen matrix
    A.resize(20,1);

    //Pointer magic
    char fileSpec[strlen(file)+1];
    snprintf( fileSpec, sizeof( fileSpec ), "%s", file );
    FILE *out = fopen( fileSpec, "r" );


    //Read and store formatted data to Matrics
    for (int i = 0; i < 20; ++i) {
        fscanf(out,"%lf\n",&A(i,0));
    }

    fclose(out);
}

// Design Matrix
void Jread(MatrixXd& A,  char* file){

    //Populate eigen matrix
    A.resize(16,20);

    //Pointer magic
    char fileSpec[strlen(file)+1];
    snprintf( fileSpec, sizeof( fileSpec ), "%s", file );
    FILE *out = fopen( fileSpec, "r" );


    //Read and store formatted data to Matrics
    for (int i = 0; i < 16; ++i) {
        fscanf(out,"%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,\n",
               &A(i,0), &A(i,1), &A(i,2), &A(i,3), &A(i,4),
               &A(i,5), &A(i,6), &A(i,7), &A(i,8), &A(i,9),
               &A(i,10), &A(i,11), &A(i,12), &A(i,13), &A(i,14),
               &A(i,15), &A(i,16), &A(i,17), &A(i,18), &A(i,19));
    }

    fclose(out);
}

// Better Print
void print(MatrixXd M, char* name, int units){
    char* u;
    switch (units) {
        case 1:
            u = "unitless";
            break;
        case 10000:
            u = "cm^2";
            break;
        default:
            cout << "Samir, you are breaking the car!" << endl;
            break;
    }
    printf("This is the %s matrix (%s)\n\n", name, u);
    int size = M.rows();
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
                if (M(i, j) > 0) {
                    printf(" %.2f ", units * M(i, j));
                } else {
                    printf("%.2f ", units * M(i, j));
                }
        }
        printf("\n");
    }
    printf("\n");
}

// Correlation Matrix
MatrixXd corr(MatrixXd C){
    int size = C.rows();
    MatrixXd R;
    R.resize(size, size);
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            if (i / (j + .5) <= 1) {
                if (i == j) {
                    R(i, j) = 1;
                } else {
                    R(i, j) = C(i, j) / (sqrt(C(i, i)) * sqrt(C(j, j)));
                    R(j, i) = R(i, j);
                }
            }
        }
    }
    return(R);
}



// Print to Files
void experimentalprint(MatrixXd M, int units, char* file){

    // Write to File
    char fileSpec[strlen(file)+1];
    snprintf( fileSpec, sizeof( fileSpec ), "%s", file );
    FILE *out = fopen( fileSpec, "w" );

    int rows = M.rows();
    int cols = M.cols();
    // Correlation/Design
    if (units == 1) {
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols-1; ++j) {
                if (abs(units * M(i,j))==1 || M(i,j)==0){
                    fprintf(out, "\\cellcolor{white}%.2f & ", units * M(i, j));
                }
                else{
                    fprintf(out, "\\cellcolor{green!%f}%.2f & ", units * M(i,j) * 100, units * M(i, j));
                }
            }

            if (abs(units * M(i,cols-1))==1 || M(i,cols-1)==0){
                fprintf(out, "\\cellcolor{white}%.2f \\\\\n\\hline\n", units * M(i, cols-1));
            }
            else{
                fprintf(out, "\\cellcolor{green!%f}%.2f \\\\\n\\hline\n",
                        units * M(i,cols-1) * 100, units * M(i, cols-1));
            }
        }
        fprintf(out,"\n");
    }
    // Covariance
    else if (units == 10000){
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                fprintf(out, "\\cellcolor{green!%f %.2f ", 10.3*units, units * M(i, j));
            }
            printf("\n");
        }
        printf("\n");
    }
    fclose(out);
}