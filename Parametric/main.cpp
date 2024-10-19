#include "header.h"

int main() {

    char* file = "C:/Users/gerge/CLionProjects/Lab_4/ctrlPnts_2024.txt";
    MatrixXd c;
    Uread(c,file);
    cout << "This is the known matrix (c) (m)" << endl << endl << c << endl << endl;


    MatrixXd x(2,1);
    x << 180, 175;
    cout << "This is the initial unknown guesses matrix (x) (m)" << endl << endl << x << endl << endl;


    file = "C:/Users/gerge/CLionProjects/Lab_4/distObs_2024.txt";
    MatrixXd l;
    Lread(l,file);
    cout << "This is the observations vector (l) (m)" << endl << endl << l << endl << endl;


    MatrixXd Cl;
    Cl.setZero(11,11);
    for (int i = 0; i < 11; ++i) {
        Cl(i,i) = pow(0.001+l(i,0)/500000, 2);
    }
    print(Cl, "Cl", 1000000);

    MatrixXd P = Cl.inverse();


    MatrixXd delta; // Corrections Vector
    MatrixXd A; // Design Matrix
    MatrixXd w; // Misclosure Vector
    MatrixXd N; // Normal Equations Matrix
    int n = 11;
    int u = 2;

    parametric(delta,0.000001,A,x,c,l,w,N,P,n,u);


    //post loop adjustment

    // Residuals Vector
    MatrixXd v = A*delta+w;
    cout << "This is the residuals vector (v) (m)" << endl << endl << v << endl << endl;

    // Write to File
    ofstream files("C:/Users/gerge/CLionProjects/Lab_4/v.txt");
    files << v;
    files.close();


    MatrixXd lhat = l+v;

    // Verify misclosure to be zero
    for (int i = 0; i < n; ++i) {
        cout << lhat(i,0) - sqrt(pow(x(0,0)-c(i,0), 2) +
                                 pow(x(1,0)-c(i,1), 2) ) << endl;
    }



    double apost = (v.transpose()*P*v).value() / (n - u);
    MatrixXd Cxhat = apost*N.inverse();
    MatrixXd Clhat = A*Cxhat*A.transpose();
    MatrixXd Cvhat = apost*Cl-Clhat;

    MatrixXd Rxhat = corr(Cxhat);
    MatrixXd Rlhat = corr(Clhat);
    MatrixXd Rvhat = corr(Cvhat);


    return 0;
}

/*
 * std for Cl
 * whats  going into
 * xhat = l naught + delta
 * non linear why?
 * initail aprox.
 * this happents this causes thiss this makes sence went interpreting correlation matrix
 * equations = measurements.
 *
 * pre adjust
 * m = n parametrix ec  quations
 * buildCl choose apriori = 1
 * build lobs
 * p = 1/preori Cl^-1
 * x nought initial guess given.
 *
 * start whil loop (delta > threashold)
 * A = df/dx nought
 * l nought = f(xnought,c)
 * w = lnought - l
 * N = At P A
 * U = At P w
 * two decimals below specified precision is good thershold
 * 
 * loop
 * Delta hat = -N^-1 u
 * abs of delta < thershold
 * xhat = xnought + delta
 *
 * less than threshold? x nought = x hat iterate again
 * good? break xhat is final
 *
 * post loop adjustment
 * vhat = A deltahat +w
 * lhat = l + vhat
 * apost = vt P v / n - u
 * Cxhat = post N^-1    this
 * Clhat = A Cxhat At   and this use for precision assurance
 * Cvhat =  Cl - Clhat  blunder detection outlier>
 *
 *
 * correcley do A and l nought rest is porcedure to follw
 * load measurements
 * loup build Cl, l, P
 * while delta is large
 * for build l nought, w, A
 * N = At P A
 * V - At P w
 * deltahat = -N^-1 U
 * xhat = x nought + delta
 * xnought = xhat
 * end
 * porst adjustmeans
 *
 */
