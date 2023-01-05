#include <string>
#include <iostream>
#include <vector>
#include <cmath>
#include <typeinfo>


#include "sparse.hpp"
#include "COO2CSR.hpp"
#include "matvecops.hpp"


/* Method to modify sparse matrix dimensions */
void SparseMatrix::Resize(int nrows, int ncols){
    j_idx.resize((nrows*ncols), 0); // Resize row index vector
    i_idx.resize((nrows*ncols), 0); // Resize col index vector
    a.resize((nrows*ncols), 0); // Resize values vector
    for(int i = 0; i < nrows; i++){
        for(int j=0; j< (ncols); j++){
            i_idx[((ncols*i) + j)] = i; // Set index for rows
            j_idx[((ncols*i) + j)] = j; // Se index for cols
        }
    }
}


void SparseMatrix::SetSize(int nrows, int ncols){
    this-> nrows = nrows; // Set value for nrows attribute
    this-> ncols = ncols; // Set value for ncols attribute
}

void SparseMatrix::AddEntry(int i, int j, double val){
    a[((ncols*i) + j)] = val; // Add value on the corresponding position
}

void SparseMatrix::ConvertToCSR(){
    COO2CSR(a, i_idx,j_idx);  
}

void SparseMatrix::PrintCOO(){
    int row;
    int col;
     std::cout << "Size a = " << a.size() << std::endl;
    for(unsigned int i = 0; i< i_idx.size(); i++){
        row = i_idx[i];
        col = j_idx[i];
        std::cout << "(" << row << "," << col << "): " << a[i] << std::endl;
    }
}

/* Method to perform sparse matrix vector multiplication using CSR formatted matrix */
std::vector<double> SparseMatrix::MulVec(std::vector<double> &vec){
    // Number of rows in the A matrix
    unsigned int n_elem = ((int)i_idx.size() -1); 
    // Resulting vector from the A*x product
    std::vector<double> y(n_elem); 
    // Initialize vector with zeros
    std::fill(y.begin(), y.end(), 0); 

    // Compute the product between each value in the CSR A matrix and the given vector.
    for (unsigned int i = 0; i < n_elem; i++){
        for (unsigned int j = (unsigned int)i_idx[i]; j < (unsigned int)i_idx[i+1]; j++){
            y[i] += a[j]*vec[j_idx[j]];
        }
    }
    return y;
}
