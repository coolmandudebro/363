#include "header.h"

// Observations Matrix
void read(MatrixXd& A,  char* file, int rows, int cols){

    //Populate eigen matrix
    A.setZero(rows,cols);

    //Pointer magic
    char fileSpec[strlen(file)+1];
    snprintf( fileSpec, sizeof( fileSpec ), "%s", file );
    FILE *out = fopen( fileSpec, "r" );

    //Read and store formatted data to Matrics
    switch (cols) {
        case 1:
            for (int i = 0; i < rows; ++i) {
                fscanf(out,"%lf\n", &A(i,0));
            }
            break;
        case 2:
            for (int i = 0; i < rows; ++i) {
                fscanf(out,"%lf %lf\n", &A(i,0), &A(i,1));
            }
            break;
    }

    fclose(out);
}

// Scientific Wild Ass Guess
void swag(MatrixXd& x, MatrixXd l){

    //find midpoint between point 1 and 25
    double x1 = ( l(0,0)+l(26,0) ) / 2;
    double y1 = ( l(0,1)+l(26,1) ) / 2;

    //find midpoint between point 1 and 75
    double x2 = ( l(0,0)+l(76,0) ) / 2;
    double y2 = ( l(0,1)+l(76,1) ) / 2;

    //find slope between point 1 and 25
    double m1 = ( l(0,1)-l(26,1) ) / ( l(0,0)-l(26,0) );

    //find slope between point 1 and 75
    double m2 = ( l(0,1)-l(76,1) ) / ( l(0,0)-l(76,0) );

    //find y-intercept between point 1 and 25
    double b1 = y1+m1*x1;

    //find y-intercept between point 1 and 75
    double b2 = y2+m2*x2;

    //find coordinates of central intercept
    x(0,0) = (b1-b2) / (m1-m2);
    x(0,1) = -m1*x(0,0)+b1;

    // radius: distance from centre to edge
    x(0,2) = sqrt( pow( x(0,0)-l(0,0) , 2) + pow( x(0,1)-l(0,1) , 2) );
}

// Better Printf
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

MatrixXd designA(int m, int u, MatrixXd l, MatrixXd x){
    MatrixXd A;
    A.setZero(m,u);
    // Populate Design Matrix
    for (int i = 0; i < m; ++i) {
        // Deviratives in terms of unknowns
        A(i,0) = -2*( l(i,0)-x(0,0) );
        A(i,1) = -2*( l(i+102,0)-x(0,1) );
        A(i,2) = -2*x(0,2);
    }
    return(A);
}

MatrixXd designB(int m, int n, MatrixXd l, MatrixXd x){
    MatrixXd B;
    B.setZero(m,n);
    // Populate Design Matrix
    for (int i = 0; i < m; ++i) {
        // Deviratives in terms of observations
        B(i,i) = 2*( l(i,0)-x(0,0) );
        B(i,i+102) = 2*( l(i+102,0)-x(0,1) );
    }
    return(B);
}

MatrixXd funcModel(MatrixXd l, MatrixXd x){
    MatrixXd dumb(102,1);
    for (int i = 0; i < 102; ++i) {
        dumb(i,0) = pow( l(i,0) - x(0,0) ,2)
                + pow( l(i+102,0) - x(0,1) ,2) - pow(x(0,2),2);
    }
    return(dumb);
}

void implicit(MatrixXd& delta, MatrixXd& v, double thresh, MatrixXd& A, MatrixXd& B, MatrixXd& x, MatrixXd l,
              MatrixXd& ln, MatrixXd& w, MatrixXd& N, MatrixXd& M, MatrixXd P, int n, int u){

    int m = n/2;

    // More Memory Allocation
    delta.resize(u,1);
    delta << 1, 1, 1;
    A.resize(m,u); // Design Matrix
    B.resize(m,n); // Design Matrix
    v.setZero(n,1); // Residuals
    w.resize(m,1); // Misclosure Vector
    N.resize(u,u); // Normal Equations Matrix
    M.resize(m,m); // M

    MatrixXd U(u,1); // Normal Equations Vector
    MatrixXd K(m,1);
    MatrixXd dumb;
    dumb = l;
    l.resize(n,1);
    l << dumb.col(0), dumb.col(1);
    MatrixXd vn(n,1);
    MatrixXd dv;
    dv.setOnes(n,1);

    int iterate = 1;

    ofstream it("C:/Users/gerge/CLionProjects/Lab_5/iterations.txt");

    // Begin
    while( delta.maxCoeff() > thresh || dv.maxCoeff() > thresh ) {
        vn = v;
        ln = l + v;

        //cout << ln << endl << endl;

        cout << iterate << endl << endl;
        it << "For iteration: " << iterate << endl << endl;

        A = designA(m,u,ln,x);
        it << "This is the design matrix (A) (unitless)" << endl << endl << A << endl << endl;

        B = designB(m,n,ln,x);
        it << "This is the design matrix (B) (unitless)" << endl << endl << B << endl << endl;

        w = funcModel(ln,x) - B*(ln-l);
        M = B*P.inverse()*B.transpose();
        N = A.transpose()*M.inverse()*A;
        U = A.transpose()*M.inverse()*w;
        delta = -N.inverse()*U;
        K = M.inverse()*(A*delta+w);
        v = -P.inverse()*B.transpose()*K;
        x = x + delta.transpose();// only transposed because i decided to make x a column vector
        dv = v - vn;
        iterate++;
        it << "This is the corrections to the unknown parameters (m)" << endl << endl << delta << endl << endl;
        it << "This is the best estimate of the unknown parameters (m)" << endl << endl << x << endl << endl;
        it << "This is the corrections to the residuals (m)" << endl << endl << v << endl << endl;
        /*cout << w << endl << endl;
        cout << delta << endl << endl;
        cout << v << endl << endl;
        cout << dv << endl << endl;*/
    }
    it.close();
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

// Latex Matrix Print
void square_matrix(MatrixXd M, char* name1, char* name2, bool hat, int units, FILE *output){

    int size = M.rows();
    if (hat){
        fprintf(output,"\\[\n%s_{\\hat{%s}\\ %d \\times %d} = \n\\]\n\\[\\begin{vmatrix*}[r]\n",
                name1, name2, size, size);
    }
    else {
        fprintf(output,"\\[\n%s_{%s\\ %d \\times %d} = \n\\]\n\\[\\begin{vmatrix*}[r]\n",
                name1, name2, size, size);
    }
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size-1; ++j) {
            fprintf(output, "%.2f & ", units*M(i,j));
        }
        fprintf(output, "%.2f \\\\ \n", units*M(i,size-1));
    }
    fprintf(output,"\\end{vmatrix*}\n\\]\n\n");

}

// Latex Table Print
void vectorf(MatrixXd M, char* name1, char* name2, bool hat, int units, FILE *output){

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

    int rows = M.rows();
    int cols = M.cols();
    // Generates a cols number of columns
    string a;
    for(int i = 0; i < cols; ++i){
        a = a + " &";
    }
    char* an = a.data();

    string t;
    for(int i = 0; i < cols; ++i){
        t = t + " r | ";
    }
    char* ti = t.data();

    if (hat){
        fprintf(output,"\\begin{tabular}{ | c |%s}\n"
                       "\\hline\n"
                       " %s %s_{\\hat{%s}_{(%s)}} \\\\ \n"
                       "\\hline\n", ti, an, name1, name2, u);
    }
    else{
        fprintf(output,"\\begin{tabular}{ | c |%s}\n"
                       "\\hline\n"
                       " %s $%s_{%s_{(%s)}}$ \\\\ \n"
                       "\\hline\n", ti, an, name1, name2, u);
    }
    for (int i = 0; i < rows; ++i) {
        fprintf(output, "%d &", i+1);
        for (int j = 0; j < cols-1; ++j) {
            fprintf(output, " %.3f &", units*M(i,j));
        }
        fprintf(output, " %.3f \\\\ \n\\hline\n", units*M(i,cols-1));
    }
    fprintf(output,"\\end{tabular}\n\n");


}
