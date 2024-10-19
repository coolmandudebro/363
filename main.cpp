#include "header.h"

int main() {

    //Write results to files
    FILE * output;
    output = fopen("/home/matyas/CLionProjects/Lab_3/results.txt","w");

    // Design Matrix import
    MatrixXd B;

    char *file = "/home/matyas/CLionProjects/Lab_3/B.csv";

    Bread(B, file);

    cout << "This is the design matrix (B) (Unitless)" << endl << endl << B << endl << endl;

    fprintf(output,"\\[\nB_{9 \\times 20} = \n\\]\n\\[\\left|\\begin{smallmatrix}\n");
    for (int i = 0; i < 9; ++i) {
        for (int j = 0; j < 20-1; ++j) {
            fprintf(output, "%.0f & ",B(i,j));
        }
        fprintf(output, "%.0f \\\\ \n",B(i,19));
    }
    fprintf(output,"\\end{smallmatrix}\\right|\n\\]\n\n");


    // Variance of Observations Matrix import
    MatrixXd Cl;

    file = "/home/matyas/CLionProjects/Lab_3/stdev_2024.txt";

    Cread(Cl, file);

    cout << "This is the variance-covariance matrix of observation (Cl) (cm^2)" << endl << endl << Cl * 10000 << endl << endl;

    fprintf(output,"\\[\nC_{l 20 \\times 20} = \n\\]\n\\[\\left|\\begin{smallmatrix}\n");
    for (int i = 0; i < 20; ++i) {
        for (int j = 0; j < 20-1; ++j) {
            if (Cl(i,j)==0){
                fprintf(output, "%.0f & ",10000*Cl(i,j));
            }
            else{
                fprintf(output, "%.2f & ",10000*Cl(i,j));
            }
        }
        if (Cl(i,19)==0) {
            fprintf(output, "%.0f \\\\ \n", Cl(i, 19));
        }
        else{
            fprintf(output, "%.2f \\\\ \n", 10000*Cl(i, 19));
        }
    }
    fprintf(output,"\\end{smallmatrix}\\right|\n\\]\n\n");


    // Observations Vector import
    MatrixXd l;

    file = "/home/matyas/CLionProjects/Lab_3/dh_2024.txt";

    Lread(l, file);

    cout << "This is the observations vector (l) (m)" << endl << endl << l << endl << endl;


    // Weight Matrix
    // A-priori variance factor = 1
    MatrixXd P = 1 * Cl.inverse();

    cout << "This is the weight matrix (cm^-2)" << endl << endl << P / 10000 << endl << endl;

    fprintf(output,"\\[\nP_{ 20 \\times 20} = \n\\]\n\\[\\left|\\begin{smallmatrix}\n");
    for (int i = 0; i < 20; ++i) {
        for (int j = 0; j < 20-1; ++j) {
            if (P(i,j)==0){
                fprintf(output, "%.0f & ",P(i,j));
            }
            else{
                fprintf(output, "%.2f & ",P(i,j)/10000);
            }
        }
        if (P(i,19)==0) {
            fprintf(output, "%.0f \\\\ \n", P(i, 19));
        }
        else{
            fprintf(output, "%.2f \\\\ \n", P(i, 19) / 10000);
        }
    }
    fprintf(output,"\\end{smallmatrix}\\right|\n\\]\n\n");

    // Misclosure Vector
    MatrixXd w = B * l;

    cout << "This is the misclosure vector (w) (mm)" << endl << endl << 1000*w << endl << endl;

    fprintf(output,"\\begin{tabular}{ | c | c | }\n"
                   "\\hline\n"
                   "Equation g_l & Misclosure (mm) \\\\ \n"
                   "\\hline\n");
    for (int j = 0; j < 9; ++j) {
        fprintf(output, "%d & %.0f \\\\ \n\\hline\n", j+1, 1000*w(j,0));
    }
    fprintf(output,"\\end{tabular}\n\n");


    // Residuals Vector
    MatrixXd v = -Cl * B.transpose() * (B * Cl * B.transpose()).inverse() * w;

    // Write to File
    ofstream files("/home/matyas/CLionProjects/Lab_3/v.txt");
    files << v;
    files.close();

    cout << "This is the residuals vector (v) (mm)" << endl << endl << 1000*v << endl << endl;

    fprintf(output,"\\begin{tabular}{ | c | c | }\n"
                   "\\hline\n"
                   "Observations \\Delta h & Residuals (mm) \\\\ \n"
                   "\\hline\n");
    for (int j = 0; j < 20; ++j) {
        fprintf(output, "%d & %.2f \\\\ \n\\hline\n", j+1, 1000*v(j,0));
    }
    fprintf(output,"\\end{tabular}\n\n");


    // A-posteriori variance factor
    double apost = (v.transpose() * P * v / 9).value();

    cout << "This is the a-posteriori variance factor: " << apost << endl << endl;;


    // Residuals Variance-Covariance Matrix
    // Aporstiriori used here
    MatrixXd Cv = apost * Cl * B.transpose() * (B * Cl * B.transpose()).inverse() * B * Cl;

    print(Cv, "Cv", 10000);

    fprintf(output,"\\[\nC_{v 20 \\times 20} = \n\\]\n\\[\\left|\\begin{smallmatrix}\n");
    for (int i = 0; i < 20; ++i) {
        for (int j = 0; j < 20-1; ++j) {
            fprintf(output, "%.2f & ",10000*Cv(i,j));
        }
        fprintf(output, "%.2f \\\\ \n",10000*Cv(i,19));
    }
    fprintf(output,"\\end{smallmatrix}\\right|\n\\]\n\n");


    // Adjusted Observations
    MatrixXd lhat = l + v;

    cout << "This is the adjusted observations vector (lhat) (m)" << endl << endl << l << endl << endl;

    fprintf(output,"\\begin{tabular}{ | c | c | }\n"
                   "\\hline\n"
                   " & Adjusted Observations \\hat{l} (m) \\\\ \n"
                   "\\hline\n");
    for (int j = 0; j < 20; ++j) {
        fprintf(output, "%d & %.3f \\\\ \n\\hline\n", j+1, lhat(j,0));
    }
    fprintf(output,"\\end{tabular}\n\n");


    // Adjusted Observations Variance-Covariance Matrix
    // Aporstiriori used here
    MatrixXd Clhat = apost*Cl - Cv;

    print(Clhat, "Clhat", 10000);

    fprintf(output,"\\[\nC_{\\hat{l} 20 \\times 20} = \n\\]\n\\[\\left|\\begin{smallmatrix}\n");
    for (int i = 0; i < 20; ++i) {
        for (int j = 0; j < 20-1; ++j) {
            fprintf(output, "%.2f & ",10000*Clhat(i,j));
        }
        fprintf(output, "%.2f \\\\ \n",10000*Clhat(i,19));
    }
    fprintf(output,"\\end{smallmatrix}\\right|\n\\]\n\n");

    // Check misclosure is numerically zero
    cout << "Check misclosure (m)" << endl << endl << B*lhat << endl << endl;


    // Correlation of Adjustment Residuals
    MatrixXd Rv = corr(Cv);

    print(Rv, "Rv", 1);

    fprintf(output,"\\[\nR_{\\hat{v} 20 \\times 20} = \n\\]\n\\[\\left|\\begin{smallmatrix}\n");
    for (int i = 0; i < 20; ++i) {
        for (int j = 0; j < 20-1; ++j) {
            fprintf(output, "%.2f & ",Rv(i,j));
        }
        fprintf(output, "%.2f \\\\ \n",Rv(i,19));
    }
    fprintf(output,"\\end{smallmatrix}\\right|\n\\]\n\n");


    // Correlation of Adjusted Observations
    MatrixXd Rlhat = corr(Clhat);

    print(Rlhat, "Rlhat", 1);

    fprintf(output,"\\[\nR_{\\hat{l} 20 \\times 20} = \n\\]\n\\[\\left|\\begin{smallmatrix}\n");
    for (int i = 0; i < 20; ++i) {
        for (int j = 0; j < 20-1; ++j) {
            fprintf(output, "%.2f & ",Rlhat(i,j));
        }
        fprintf(output, "%.2f \\\\ \n",Rlhat(i,19));
    }
    fprintf(output,"\\end{smallmatrix}\\right|\n\\]\n\n");




    /////////// Task 2 /////////////////



    // Import Jacobian Matrix
    MatrixXd J;

    file = "/home/matyas/CLionProjects/Lab_3/J.csv";

    Jread(J, file);

    cout << "This is the jacobian matrix (J) (unitless)" << endl << endl << J << endl << endl;

    fprintf(output,"\\[\nJ_{16 \\times 20} = \n\\]\n\\[\\left|\\begin{smallmatrix}\n");
    for (int i = 0; i < 16; ++i) {
        for (int j = 0; j < 20-1; ++j) {
            fprintf(output, "%.0f & ",J(i,j));
        }
        fprintf(output, "%.0f \\\\ \n",J(i,19));
    }
    fprintf(output,"\\end{smallmatrix}\\right|\n\\]\n\n");

    
    // Preload Base station height
    MatrixXd height(16,1);
    height.fill(148.467);
    
    // Populate Unknown Parameters Matrix
    MatrixXd x = height+(J*lhat);

    cout << "This is the unknown parameters vector (x) (m)" << endl << endl << x << endl << endl;

    fprintf(output,"\\begin{tabular}{ | c | c | }\n"
                   "\\hline\n"
                   " & Unknown Parameters \\hat{x} (m) \\\\ \n"
                   "\\hline\n");
    for (int j = 0; j < 16; ++j) {
        fprintf(output, "%d & %.3f \\\\ \n\\hline\n", j+1, x(j,0));
    }
    fprintf(output,"\\end{tabular}\n\n");


    // Unknown Parameters Variance-Covariance Matrix
    MatrixXd Cx = J*Clhat*J.transpose();

    print(Cx, "Cx",10000);

    fprintf(output,"\\[\nC_{\\hat{x} 16 \\times 16} = \n\\]\n\\[\\left|\\begin{smallmatrix}\n");
    for (int i = 0; i < 16; ++i) {
        for (int j = 0; j < 16-1; ++j) {
            fprintf(output, "%.2f & ",10000*Cx(i,j));
        }
        fprintf(output, "%.2f \\\\ \n",10000*Cx(i,15));
    }
    fprintf(output,"\\end{smallmatrix}\\right|\n\\]\n\n");


    // Correlation of Unknown Parameters
    MatrixXd Rxhat = corr(Cx);

    print(Rxhat, "Rxhat",1);

    fprintf(output,"\\[\nR_{\\hat{x} 20 \\times 20} = \n\\]\n\\[\\left|\\begin{smallmatrix}\n");
    for (int i = 0; i < 16; ++i) {
        for (int j = 0; j < 16-1; ++j) {
            fprintf(output, "%.2f & ",Rxhat(i,j));
        }
        fprintf(output, "%.2f \\\\ \n",Rxhat(i,15));
    }
    fprintf(output,"\\end{smallmatrix}\\right|\n\\]\n\n");

    fclose(output);

    return (0);
}



/*
 methodology is theory not i did this that
 do paper
 do text file
 aposteriora variance factor use as scalf factor for variance covariance corraelation matices
 residuals tell you precission
 what is small residual? what is your presission; specs, relative?
 part 2 seperate from adjustment
 condition equation: f = delta h(2-1) + delta h(7-2) - delta h(7-1) = 0
 linearly indipendant paths
 connected linearly independant
 visualize design matrix
 rows f1, f2, f3, .....
 columns h(21) h 22, h23, h24
 tabs to seperate b matrix
 check very
 look at direction of arrow for
 Cl = NxN = variance and zeros
 use 0 as initial estimate for easy checking
 design weight misclosur
 rank is less than redundancey therefore dependant, fix
 matlab rank function
 lhat = l + vhat
 h2 = delta h b + delta h
 build jacobian = v x l =
 cx = estimated precision of unknown
 redundancey = measuarements - unknowns. higher is better to check if ans is good
 scale with aposterior factor before claculate Cvhat
 adjusted measurements could be correlated
 aposteriora variance factor simga hat naught ^2 = Vt p V / r if >1 too optemistic, if <1 too pecimistic, 1 = good, maybe depends
 read what the lab is asking for referance all outputs
*/