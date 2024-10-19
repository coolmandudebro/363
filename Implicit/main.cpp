#include "header.h"

int main() {

    char* file = "C:/Users/gerge/CLionProjects/Lab_5/XY_coords_2024.txt";
    char* out = "C:/Users/gerge/CLionProjects/Lab_5/output.txt";
    FILE *output = fopen( out, "w" );
    ofstream results("C:/Users/gerge/CLionProjects/Lab_5/results.txt");


    // Observations matrix
    MatrixXd l;
    read(l,file,102,2);
    //cout << "This is the observations vector (l) (m)" << endl << endl << l << endl << endl;
    results << "This is the observations vector (l) (m)" << endl << endl << l << endl << endl;
    //vectorf(l,"l","", false,1,output);


    // Unknowns vector
    MatrixXd x(1,3);
    swag(x,l);
    cout << "This is the initial unknown guesses matrix (x) (m)" << endl << endl << x << endl << endl;
    results << "This is the initial unknown guesses matrix (x) (m)" << endl << endl << x << endl << endl;
    vectorf(x,"x","", false,1,output);


    // Variance of the Observations
    MatrixXd Cl;
    Cl.setZero(204,204);
    for (int i = 0; i < 204; ++i) {
        Cl(i,i) = 0.000009;
    }
    results << "This is the Variance of the Observations matrix (Cl) (mm^2)" << endl << endl << 1000000*Cl << endl << endl;
    /*print(Cl, "Cl", 1000000);
    square_matrix(Cl,"C","l",false,1000000,output);*/


    // A-priori variance factor of 1.
    MatrixXd P = 1*Cl.inverse();
    results << "This is the weight matrix (P) (mm^-2)" << endl << endl << P/1000000 << endl << endl;


    // Initalizing Memory for the Adjustment
    MatrixXd delta; // Corrections Vector
    MatrixXd v;
    MatrixXd A; // Design Matrix
    MatrixXd B; // Design Matrix
    MatrixXd w; // Misclosure Vector
    MatrixXd N; // Normal Equations Matrix
    MatrixXd M; // M
    MatrixXd lhat; // Adjusted Observations
    int n = 204; // Observations
    int u = 3; // Unknowns


    // Adjustment Iterative Loop
    implicit(delta,v,0.00001,A,B,x,l,lhat,w,N,M,P,n,u);


    // // // // // // // Post loop adjustment // // // // // // //


    // Write to File for matlab plot
    ofstream files("C:/Users/gerge/CLionProjects/Lab_5/v.txt");
    files << v;
    files.close();

    //cout << "This is the vector of residuals (vhat) (mm)" << endl << endl << 1000*v << endl << endl;
    results << "This is the vector of residuals (vhat) (mm)" << endl << endl << 1000*v << endl << endl;

    cout << "This is the vector of best estimates (xhat) (m)" << endl << endl << x << endl << endl;
    results << "This is the vector of best estimates (xhat) (m)" << endl << endl << x << endl << endl;
    vectorf(x,"x","", true,1,output);

    //cout << "This is the vector of adjusted observations (lhat) (m)" << endl << endl << lhat << endl << endl;
    results << "This is the vector of adjusted observations (lhat) (m)" << endl << endl << lhat << endl << endl;

    // Aposteriori Variance Factor
    double apost = (v.transpose()*P*v).value() / (n - u);
    cout << "This is the a-posteriori variance factor: " << apost << endl << endl;
    results << "This is the a-posteriori variance factor: " << apost << endl << endl;

    // Variance-Covariance of the Unknowns
    MatrixXd Cxhat = apost*N.inverse();
    print(Cxhat, "Cxhat", 1000000);
    results << "This is the Variance-Covariance matrix of the Unknowns (Cxhat) (mm^2)" << endl << endl << 1000000*Cxhat << endl << endl;
    square_matrix(Cxhat,"C","x",true,1000000,output);

    
    // Variance-Covariance of the Residuals
    MatrixXd Cvhat = apost*P.inverse()*B.transpose()*M.inverse()*B*P.inverse()
            - apost*P.inverse()*B.transpose()*M.inverse()*A*N.inverse()*A.transpose()*M.inverse()*B*P.inverse();
    results << "This is the Variance-Covariance matrix of the Residuals (Cvhat) (mm^2)" << endl << endl << 1000000*Cvhat << endl << endl;
    /*print(Cvhat, "Cvhat", 1000000);
    square_matrix(Cvhat,"C","v",true,1000000,output);*/


    // Variance-Covariance of the Adjusted Observations
    MatrixXd Clhat = apost*Cl-Cvhat;
    results << "This is the Variance-Covariance matrix of the Adjusted Observations (Clhat) (mm^2)" << endl << endl << 1000000*Clhat << endl << endl;
    /*print(Clhat, "Clhat", 1000000);
    square_matrix(Clhat,"C","l",true,1000000,output);*/


    // Check equation fulfillment
    cout << "This is the misclosure of the functional model (m)" << endl << endl << funcModel(lhat,x) << endl << endl;
    results << "This is the misclosure of the functional model (m)" << endl << endl << funcModel(lhat,x) << endl << endl;

    // Correlation Coefficients of the Unknowns
    MatrixXd Rxhat = corr(Cxhat);
    print(Rxhat, "Rxhat", 1);
    square_matrix(Rxhat,"R","x",true,1,output);
    results << "This is the Correlation Coefficient matrix of the Unknowns (Rxhat) (unitless)" << endl << endl << Rxhat << endl << endl;

    // Correlation Coefficients of the Adjusted Observaions
    MatrixXd Rlhat = corr(Clhat);
    /*print(Rlhat, "Rlhat", 1);
    square_matrix(Rlhat,"R","l",true,1,output);*/
    results << "This is the Correlation Coefficient matrix of the Adjusted Observations (Rlhat) (unitless)" << endl << endl << Rlhat << endl << endl;

    // Correlation Coefficients of the Residuals
    MatrixXd Rvhat = corr(Cvhat);
    /*print(Rvhat, "Rvhat", 1);
    square_matrix(Rvhat,"R","v",true,1,output);*/
    results << "This is the Correlation Coefficient matrix of the Residuals (Rvhat) (unitless)" << endl << endl << Rvhat << endl << endl;


    // Close files
    fclose(output);
    results.close();

    return 0;
}
/*
 * subplot function matlab one in each quadrent
 *
 * Combined model loop.
 *
 * Pre Adjustment
 *
 * m = n/2sb = 102 -> f(x,l,c) = 0
 * r = m - u = 99
 * aprior = 1 happy //happy? wym?
 * Build Cl, lobs (diff from ln)
 * xn = ? (u x 1)
 * ln = (not initial values, not parametric) = l
 * A = d f / d x, at xn and ln (m x u)
 * B = d f / d l, at xn and ln (m x n)
 *
 * end pre adjustment
 *
 * Begin Loop
 * while ( |delta (d) , v| > thresh
 * Build
 * B = d f / d ln
 * A = d f / d xn
 * w = (is how did the equation close) = f(xn,ln,c) - B(ln-l)
 * jacobian linear transformation from the equation space to the misclosure space
 * M = B P^-1 B^T
 * N = A^T M^-1 A
 * U = A^t M^-1 w
 * dhat = -N^-1 U
 * khat = M^-1 ( A dhat+ w )
 * vhat = -P^-1 B^T k
 * xhat = xn + dhat
 * lhat = l + vhat
 * d vhat = vhat_i - vhat_{i-1}
 * end of loop
 * check | dhat and d vhat | < thershold
 * repeat with xn = xhat and ln = lhat
 * do again until leave
 *
 * End Loop
 *
 * Post Adjustment
 *
 * apost = normal
 * Cxhat = apost N^-1
 * Cvhat = apost P^-1 B^T M^-1 B P^-1  -  apost P^-1 B^T M^-1 A N^-1 A^T M^-1 B P^-1
 * AAHAHAAHAHAAHAHAAHAHAAHAAHAHAHAAHAH!!!!!!!!!!!!!!!!!!!!!!!
 *
 *
 * average  ofcoordinates is likely the center of cordinatets
 * or any three chrods
 *
 * range for radius
 *
 * Biomedical photogrammetry
 * Nose lab
*/