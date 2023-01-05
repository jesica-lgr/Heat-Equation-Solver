#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <typeinfo>
#include <vector>

#include "sparse.hpp"
#include "heat.hpp"
#include "CGSolver.hpp"

/* Setup Ax = b system */
int HeatEquation2D::Setup(std::string inputfile){
    unsigned int nrows_input = 0;
    unsigned int ncols_input = 0;
    double Tc;
    double Th;
    double length;
    double width;
    double h;
    int n_eqs;
    std::vector<std::tuple<unsigned int, unsigned int>> v;

    std::ifstream input_data;

    input_data.open(inputfile.c_str()); // Open input file
    if (input_data.is_open()){
        input_data >> length >> width >> h; // Get geometric parameters
        // Compute number of rows and columns
        nrows_input= (unsigned int)(width/h); 
        ncols_input = (unsigned int)(length/h);
        input_data >> Tc >> Th; // Get lower and upper boundary temperature
        input_data.close();
    } else {
        std::cerr << "Failed to open input file" << std::endl;
    }
    // Map between grid point ij and its corresponding position in x solution vector
    gridpoint_to_x = Create_GridVecMapping(x, v, n_eqs, ncols_input, nrows_input);
    A.SetSize(n_eqs, n_eqs); 
    // Fill matrix A in COO format and vector b 
    Fill_Ab(v, n_eqs, ncols_input, b, Th, Tc, length, h); 
    A.ConvertToCSR(); 
    return 0;
}

std::map<std::tuple<int, int>, int> HeatEquation2D::Create_GridVecMapping(std::vector<double> &x, std::vector<std::tuple<unsigned int, unsigned int>> &v, int &n_eqs, const int &ncols_input, const int &nrows_input){

    int order_elem = 0; // Order of variables in the solution vector.
    for(unsigned int j = 0; j < (unsigned int)(ncols_input); j++){
        if(j == (unsigned int)ncols_input -1){ // Skip right boundary values
            continue;
        }
        for(unsigned int i = 0; i< (unsigned int)nrows_input; i++){
            if((i == 0) || (i == (unsigned int)(nrows_input - 1))){ // Skip the top and buttom boundaries
                continue;
            } else { 
                // Add unknown to vector x 
                x.push_back((double)1);
                // Add the position indexes of that point in the grid (input geometry) 
                v.push_back(std::make_tuple(i,j)); 
                // Map from grid point ij to position in x
                gridpoint_to_x.insert({std::make_tuple(i,j), order_elem}); 
                order_elem = order_elem +1; // next element number
            }

        }
    }
    n_eqs = (int)x.size(); // Total number of equations
    return gridpoint_to_x;
}



void HeatEquation2D::Fill_Ab(std::vector<std::tuple<unsigned int, unsigned int>> &v, unsigned const int &n_eqs, const long unsigned int &ncols_input, std::vector<double> &b, const double &Th, const double &Tc, const double &length, const double &h){

    b.resize(n_eqs); // Set size of constant vector b
    fill(b.begin(), b.end(), 0); // Initialize b vector 

    A.Resize(n_eqs, n_eqs); // 

    std::tuple< int,  int> top;
    std::tuple< int,  int> down;
    std::tuple< int,  int> right;
    std::tuple< int,  int> left;

    for(unsigned int i = 0; i < n_eqs; i++){
        int grid_i = std::get<0>(v[i]);
        int grid_j = std::get<1>(v[i]);

        // Get the position of the point at the top and buttom of each coordinate associated with each entry in the solution's vector.
        top = std::make_tuple(grid_i-1, grid_j); 
        down = std::make_tuple(grid_i+1, grid_j);

        // Left boundary i.e. grab value from the opposite side
        if((grid_j -1) < 0){ 
            left = std::make_tuple(grid_i, ncols_input - 2);
        } else { // Consider point at the left side
            left = std::make_tuple(grid_i, grid_j - 1);
        }
        
        // Right boundary i.e. do not consider as a variable (make it equal to the first column of the grid)
        if((unsigned int)(grid_j + 1) == (unsigned int)(ncols_input - 1)){ 
            right = std::make_tuple(grid_i,  0 );
        } else { // Consider point at the right side
            right = std::make_tuple(grid_i, grid_j + 1);
        }

        // Check if point at the given direction is in the solution's vector x.
        if(gridpoint_to_x.find(top) != gridpoint_to_x.end()){ 
            A.AddEntry(i, gridpoint_to_x[top], 1); 
        } else{ // Top point is in the upper boundary
            b[i] -= Th; // Set hot temperature
        }
        if(gridpoint_to_x.find(down) != gridpoint_to_x.end()){
            A.AddEntry(i, gridpoint_to_x[down], 1); 
        } else{ // Down point is in the lower boundary
            b[i] -= Compute_Temp(grid_j, Tc, length, h); // Set cold temperature
        }
        if(gridpoint_to_x.find(right) != gridpoint_to_x.end()){
            A.AddEntry(i, gridpoint_to_x[right], 1); 
        }
        if(gridpoint_to_x.find(left) != gridpoint_to_x.end()){
            A.AddEntry(i, gridpoint_to_x[left], 1); 
        }

        A.AddEntry(i,i,-4); // Add coeficient of current unknown (i=j)
    }

}

double HeatEquation2D::Compute_Temp(const int grid_j, const double Tc, const double length, const double h){
    // Cool jet temperature function with respect to the column position
    double T = ((-1)*Tc) * (exp(-10 * pow(((grid_j*h)-(length/2)),2)) -2);
    return T;
}


int HeatEquation2D::Solve(std::string soln_prefix){
    int n_iter_status; 
    double tol = 1.e-5; // Tolerance value for solver
    // Call the Conjugate Gradient Solver
    n_iter_status = CGSolver(A.a, A.i_idx, A.j_idx, b, x,tol, soln_prefix);

    if (n_iter_status != -1){ // Solver converged. n_iter_status = number of iterations
        std::cout << "SUCCESS: CG solver converged in " << n_iter_status << " iterations." << std::endl;
        return 0; // Status evaluates to false i.e. no error message is displayed
    } else{ // Solver did not converged. n_iter_status = -1 (evaluates to True)
        return n_iter_status;
    }
}