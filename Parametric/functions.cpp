#include "header.h"

// Observations Matrix
void Lread(MatrixXd& A,  char* file){

    //Populate eigen matrix
    A.setZero(11,1);

    //Pointer magic
    char fileSpec[strlen(file)+1];
    snprintf( fileSpec, sizeof( fileSpec ), "%s", file );
    FILE *out = fopen( fileSpec, "r" );


    //Read and store formatted data to Matrics
    for (int i = 0; i < 11; ++i) {
        fscanf(out,"%lf\n", &A(i,0));
    }

    fclose(out);
}

// Unknown Matrix
void Uread(MatrixXd& A,  char* file){

    //Populate eigen matrix
    A.resize(11,2);

    //Pointer magic
    char fileSpec[strlen(file)+1];
    snprintf( fileSpec, sizeof( fileSpec ), "%s", file );
    FILE *out = fopen( fileSpec, "r" );


    //Read and store formatted data to Matrics
    for (int i = 0; i < 11; ++i) {
        fscanf(out,"%lf %lf\n", &A(i,0), &A(i,1));
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
        case 1000000:
            u = "mm^2";
            break;
        default:
            cout << "Samir, you are breaking the car!" << endl;
            break;
    }
    printf("This is the %s matrix (%s)\n\n", name, u);
    int size = M.rows();
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            if (M(i, j) >= 0) {
                printf(" %.2f ", units * M(i, j));
            } else {
                printf("%.2f ", units * M(i, j));
            }
        }
        printf("\n");
    }
    printf("\n");
}

// Adjustment Loop
void parametric(MatrixXd& delta, double thresh, MatrixXd& A, MatrixXd& x, MatrixXd c,
                MatrixXd l, MatrixXd& w, MatrixXd& N, MatrixXd P, int n, int u){
    delta.resize(u,1);
    delta << 1, 1;
    A.resize(n,u); // Design Matrix
    w.resize(n,1); // Misclosure Vector
    N.resize(u,u); // Normal Equations Matrix


    MatrixXd ln(n,1); // Approximate Observations Vector
    MatrixXd U(u,1); // Normal Equations Vector
    int iterate = 1;

    while( abs( delta(0,0) ) > thresh || abs( delta(1,0) ) > thresh ){
        cout << "for iteration: " << iterate << endl << endl;


        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < u; ++j) {
                A(i,j) = ( x(j,0)-c(i,j) ) /
                         sqrt(pow(x(j,0)-c(i,j), 2) +
                              pow(x((u-1)-j,0)-c(i,(u-1)-j), 2) );
            }
        }
        cout << "This is the design matrix (A) (unitless)" << endl << endl << A << endl << endl;


        for (int i = 0; i < n; ++i) {
            ln(i,0) = sqrt(pow(x(0,0)-c(i,0), 2) +
                           pow(x(1,0)-c(i,1), 2) );
        }
        cout << "This is the provisional observations vector (ln) (m)" << endl << endl << ln << endl << endl;


        w = ln - l;
        cout << "This is the misclosure vector (w) (m)" << endl << endl << w << endl << endl;


        N = A.transpose()*P*A;
        U = A.transpose()*P*w;
        delta = -N.inverse()*U;
        cout << "This is the corrections vector (delta) (m)" << endl << endl << delta << endl << endl;

        x = x + delta;
        cout << "This is the next unknowns vector (x) (m)" << endl << endl << x << endl << endl;
        iterate++;
    }
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

void square_matrix(MatrixXd M, char* name1, char* name2, bool hat, int units, char* file){
    //Pointer magic
    char fileSpec[strlen(file)+1];
    snprintf( fileSpec, sizeof( fileSpec ), "%s", file );
    FILE *output = fopen( fileSpec, "w" );

    int size = M.rows();
    if (hat == 1){
        fprintf(output,"\\[\n%s_{\\hat{%s} %d \\times %d} = \n\\]\n\\[\\left|\\begin{smallmatrix*}[r]\n",
                name1, name2, size, size);
    }
    else {
        fprintf(output,"\\[\n%s_{%s %d \\times %d} = \n\\]\n\\[\\left|\\begin{smallmatrix*}[r]\n",
                name1, name2, size, size);
    }
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size-1; ++j) {
            fprintf(output, "%.2f & ", units*M(i,j));
        }
        fprintf(output, "%.2f \\\\ \n", units*M(i,size-1));
    }
    fprintf(output,"\\end{smallmatrix*}\\right|\n\\]\n\n");

    fclose(output);
}

void vector(MatrixXd M, char* name1, char* name2, int units, char* file){
    //Pointer magic
    char fileSpec[strlen(file)+1];
    snprintf( fileSpec, sizeof( fileSpec ), "%s", file );
    FILE *output = fopen( fileSpec, "w" );

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

    int size = M.rows();
    fprintf(output,"\\begin{tabular}{ | c | r | }\n"
                   "\\hline\n"
                   " & %s \\hat{%s} (%s) \\\\ \n"
                   "\\hline\n", name1, name2,u);
    for (int j = 0; j < size; ++j) {
        fprintf(output, "%d & %.3f \\\\ \n\\hline\n", j+1, units*M(j,0));
    }
    fprintf(output,"\\end{tabular}\n\n");

    fclose(output);
}

