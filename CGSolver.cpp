#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <sstream>

#include "CGSolver.hpp"
#include "matvecops.hpp"



/* Solve Ax = b using the Conjugate Gradient method. */
int CGSolver(std::vector<double> &val,
             std::vector<int>    &row_ptr,
             std::vector<int>    &col_idx,
             std::vector<double> &b,
             std::vector<double> &x,
             double              tol, std::string soln_prefix){
    

    // Declare and initialize the vectors
    std::vector<double> x_0 = x; // Initial solution vector
    std::vector<double> r_0; // Stores difference of b and y
    std::vector<double> y; // A*x = y
    unsigned int n_elem = (int)x.size();
    std::vector<double> r_n(n_elem);
    std::vector<double> p_n(n_elem);
    std::vector<double> Ap_n;
    std::vector<double> x_n;
    std::vector<double> col_vec_r;
    std::vector<double> denominator_alpha;
    std::vector<double> denominator_beta;
    std::vector<double> numerator_beta;
    std::vector<double> numerator_alpha;

    y = Multiply_CSRVec(val,row_ptr,col_idx,x_0); // A*x_0 = y
    r_0 = Substract_vec(b,y); // r_0 = b -y
    std::vector<double> p_0 = r_0;

    // Compute the L2norm of initial vector
    double L2normr0 = L2norm(r_0); 
    double L2normr; 
    // Coeficients for the Conjugate Gradient Algorithm
    double beta; 
    double alpha;

    // Initial and maximum number of iterations
    unsigned int n_iter = 0;
    unsigned int n_itermax = (int)x.size();

    // While maximum number of iterations is not reached
    while (n_iter < n_itermax){
        // Transporse variables for r and p vectors
        std::vector<double> valT_r_n;
        std::vector<int> i_idxT_r_n;
        std::vector<int> j_idxT_r_n;
        std::vector<double> valT_p_n;
        std::vector<int> i_idxT_p_n;
        std::vector<int> j_idxT_p_n;

        // First iteration
        if (n_iter == 0){
            r_n = r_0;
            p_n = p_0;
            x_n = x_0;
        }

        if (n_iter%10 == 0){ // Store x every 10 iterations
            std::ofstream solution;
            std::string solution_filename;
            if (n_iter == 0){ // Initial iteration
                solution_filename = soln_prefix + "00" + std::to_string(n_iter) + ".txt";
            } else if (n_iter < 100){
                solution_filename = soln_prefix + "0" + std::to_string(n_iter) + ".txt"; 
            } else{
                solution_filename = soln_prefix + std::to_string(n_iter) + ".txt"; 
            }
            // Write the solution vector on the output file 
            solution.open(solution_filename.c_str());
            if (solution.is_open()) {
                std::ostringstream value_x;

                for (unsigned int i = 0; i < x_n.size(); i++){
                    value_x.str("");
                    value_x << std::setprecision(4) << x_n[i]; 
                    solution << value_x.str() << std::endl;
                }
            } else{
                    std::cerr << "Failed to open solutions file" << std::endl;
            }
            solution.close();
        }
        // Update number of iteration
        n_iter = n_iter + 1;
        // Compute transpose of r_n vector (r_nT).
        TransposeVec_CSR(r_n, valT_r_n, i_idxT_r_n, j_idxT_r_n); 
        numerator_alpha = Multiply_CSRVec(valT_r_n,i_idxT_r_n,j_idxT_r_n,r_n);
        // Compute transpose of p_n vector (p_nT).
        TransposeVec_CSR(p_n, valT_p_n, i_idxT_p_n, j_idxT_p_n); 
        // val, row_ptr, col_idx now define the Ap_n vector in CSR form. 
        Ap_n = Multiply_CSRVec(val,row_ptr,col_idx,p_n); 
        denominator_alpha = Multiply_CSRVec(valT_p_n,i_idxT_p_n,j_idxT_p_n,Ap_n);
        alpha = numerator_alpha[0]/denominator_alpha[0]; // Compute alpha
        std::vector<double> alpha_pn = Product_NumVec(alpha, p_n);
        // Update x_n
        x_n = Add_vec(x_n, alpha_pn); 
        std::vector<double> alpha_Apn = Product_NumVec(alpha, Ap_n);
        // Update r_n vector
        r_n = Substract_vec(r_n, alpha_Apn);
         // Compute new L2norm of r_n vector 
        L2normr = L2norm(r_n);

        // Convergence is reached. Get out of the while loop.
        if ((L2normr/L2normr0) < tol){ 
            break;
        }   

        // Convergence is not reached, update r_n parameters in CSR form for next iteration.
        denominator_beta = numerator_alpha;
        std::vector<double> valT_r_n1;
        std::vector<int> i_idxT_r_n1;
        std::vector<int> j_idxT_r_n1;


        TransposeVec_CSR(r_n, valT_r_n1, i_idxT_r_n1, j_idxT_r_n1);

        numerator_beta = Multiply_CSRVec(valT_r_n1,i_idxT_r_n1,j_idxT_r_n1,r_n);
        // Compute beta
        beta = numerator_beta[0]/denominator_beta[0]; 
        std::vector<double> beta_pn = Product_NumVec(beta, p_n);
        // Update p_n
        p_n = Add_vec(r_n,beta_pn); 

    }

    if (n_iter > n_itermax){
        return -1; // Solver did not converge
    } else {
        return (int)n_iter; // Solver converged
    }
    
}
