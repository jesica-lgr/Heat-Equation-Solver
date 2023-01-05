#ifndef MATVECOPS_HPP
#define MATVECOPS_HPP

#include <vector>
#include <iostream>

/* matvecops contains vector and matrixes common operations used in the CG algorithm. */

/*Multipy a CSR matrix A and a vector x. Returns a 'y' vector of the same size as 'x' (A*x = y). */
std::vector<double> Multiply_CSRVec(std::vector<double> &val,
             std::vector<int>    &row_ptr,
             std::vector<int>    &col_idx,
             std::vector<double> &x);


/* Substract two vectors. Returns a vector of the same size.*/
std::vector<double> Substract_vec(std::vector<double>    &vec1, 
                    std::vector<double>    &vec2);
/* Add two vectors. Returns a vector of the same size.*/                  
std::vector<double> Add_vec(std::vector<double>    &vec1,
                    std::vector<double>    &vec2);
/* Compute the product of an scalar times a vector. Returns a vector of the same size as the initial one. */
std::vector<double> Product_NumVec(double    &num,
                    std::vector<double>    &vec);

/* Computes the L2norm of a vector. Returns a double.*/
double L2norm(std::vector<double> &vec);

/* Print a vector. */
void Print_vec(std::vector<double> &vec);
void Print_vec(std::vector<int> &vec);

/* Transpose a vector and modify the CSR parameters of the transposed vector. */
void TransposeVec_CSR(std::vector<double> &vec, std::vector<double> &valT, std::vector<int> &i_idxT, std::vector<int> &j_idxT);


#endif /* MATVECOPS_HPP */