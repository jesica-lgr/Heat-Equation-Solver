#ifndef SPARSE_HPP
#define SPARSE_HPP

#include <vector>
#include <iostream>
#include <map>

#pragma once 


/* The SparseMatrix class defines a sparse matrix object, which data attributes corespond to the vectors defining the matrix. The methods contained in this class perform the necessary tasks: setsize, resize, convert COO to CSR format, add entries, print matrix, and compute the product of a matrix times a vector. */

class SparseMatrix
{
  private:
    std::vector<int> i_idx; // Row pointer vector of matrix A in CSR format
    std::vector<int> j_idx; // Col vector of matrix A in CSR format
    std::vector<double> a; // Values vector of matrix A
    int ncols; // Number of rows of A
    int nrows; // Number of columns of A
    
    friend class HeatEquation2D; 
    
  public:

    /* Method to modify sparse matrix dimensions (i_idx, j_idx, a).*/
    void Resize(int nrows, int ncols);

    /* Method to set the value for the number of rows and columns of the sparse matrix. */
    void SetSize(int nrows, int ncols);

    /* Method to add entry to matrix in COO format */
    void AddEntry(int i, int j, double val);

    /* Method to convert COO matrix to CSR format using provided function */
    void ConvertToCSR();

    /* Method to display/print a matrix in COO format. */
    void PrintCOO();

    /* Method to perform sparse matrix vector multiplication using CSR formatted matrix */
    std::vector<double> MulVec(std::vector<double> &vec);

};

#endif /* SPARSE_HPP */
