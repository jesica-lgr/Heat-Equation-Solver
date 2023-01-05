#include <vector>
#include <cmath>
#include <iostream>


#include "matvecops.hpp"
#include "COO2CSR.hpp"


std::vector<double> Multiply_CSRVec(std::vector<double> &val,
             std::vector<int>    &row_ptr,
             std::vector<int>    &col_idx,
             std::vector<double> &x){

                unsigned int n_elem = ((int)row_ptr.size() -1); // Number of rows in the A matrix
                std::vector<double> y(n_elem); // Resulting vector from the A*x product
                std::fill(y.begin(), y.end(), 0); // Initialize vector with zeros

                // Compute the product between each value in the CSR A matrix and the given vector.
                for (unsigned int i = 0; i < n_elem; i++){
                    for (unsigned int j = (unsigned int)row_ptr[i]; j < (unsigned int)row_ptr[i+1]; j++){
                        y[i] += val[j]*x[col_idx[j]];
                    }
                }
            
                return y;

}

std::vector<double> Substract_vec(std::vector<double>    &vec1,
                    std::vector<double>    &vec2){

                    unsigned int n_elem = (int)vec1.size(); 
                    std::vector<double> resultant_vec;
                    double result;

                    // Substract each entry from both vectors. Element wise substraction.
                    for (unsigned int i = 0; i < n_elem; i++){
                        result = vec1[i] - vec2[i];
                        resultant_vec.push_back(result);
                    }

                    return resultant_vec;

}

std::vector<double> Add_vec(std::vector<double>    &vec1,
                    std::vector<double>    &vec2){

                    unsigned int n_elem = (int)vec1.size();
                    std::vector<double> resultant_vec;
                    double result;

                    // Add each entry from both vectors. Element wise substraction.
                    for (unsigned int i = 0; i < n_elem; i++){
                        result = vec1[i] + vec2[i];
                        resultant_vec.push_back(result);
                    }

                    return resultant_vec;

}

std::vector<double> Product_NumVec(double    &num,
                    std::vector<double>    &vec){

                    unsigned int n_elem = (int)vec.size();
                    std::vector<double> resultant_vec;
                    double result;
                    // Element wise product of an scalar times the given vector.
                    for (unsigned int i = 0; i < n_elem; i++){
                        result = num*vec[i];
                        // Add product result to the resultant vector.
                        resultant_vec.push_back(result); 
                    }

                    return resultant_vec;

}



double L2norm(std::vector<double> &vec){
    unsigned int n_elem = (int)vec.size();
    double sum_squares = 0;
    double l2norm_result;
    
    // Square and sum each entry of the vector
    for (unsigned int i = 0; i < n_elem; i++){
        sum_squares += pow(vec[i],2);
    }
    // Compute the square root of the sum
    l2norm_result = sqrt(sum_squares);

    return l2norm_result;

}

void Print_vec(std::vector<double> &vec){
    unsigned int n_elem = (int)vec.size();

    // Output each entry of the given vector.
    for (unsigned int i = 0; i < n_elem; i++){
        std::cout << vec[i] << ", ";
    }
    std::cout << std::endl;
}

void Print_vec(std::vector<int> &vec){
    unsigned int n_elem = (int)vec.size();

    // Output each entry of the given vector.
    for (unsigned int i = 0; i < n_elem; i++){
        std::cout << vec[i] << ", ";
    }
    std::cout << std::endl;
}

void TransposeVec_CSR(std::vector<double> &vec, std::vector<double> &valT, std::vector<int> &i_idxT, std::vector<int> &j_idxT){

                std::vector<double> val;
                std::vector<int>    i_idx;
                std::vector<int>    j_idx;
                unsigned int n_elem = (unsigned int)vec.size();
        
                // Get parameters of vector 'vec' of size (1 x n_elem) in CSR form
                unsigned int non_zero_values = 0; // Stores the number of non zero values of the row

                // Inspect each element in the given vector.
                for (unsigned int j = 0; j < n_elem; j++){
                    if (vec[j] == 0){ // Column entry is a zero
                        continue;
                    }
                    non_zero_values += 1; // Count column entry value different than zero
                    valT.push_back(vec[j]); // Add column entry value to the values parameter
                    //i_idxT.push_back(j+1);
                    j_idxT.push_back(j); // Add column index to the columns parameter
                }

                // For the transpose vector, row pointer has just two values (number of element at the beginning = 0 and number of elements at the end of the (1 x n_elem) vector. 
                i_idxT.push_back(0); 
                i_idxT.push_back(non_zero_values);
                

    }
