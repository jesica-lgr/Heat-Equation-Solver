#ifndef HEAT_HPP
#define HEAT_HPP

#include "sparse.hpp"

/* HeatEquation2D contains the data and methods required to from the system of equations from a given geometry (filename), solve the system using the CGSolver and output the result to files.*/

class HeatEquation2D
{
  private:
    SparseMatrix A; // Sparse matrix
    std::vector<double> b, x; // Solution vector (x) and constants (b)
    int n_eqs; // Number of equations
    // Map between the grid/geometry coordinates and position in x vector
    std::map<std::tuple<int, int>, int> gridpoint_to_x; 
    // Vector containing the grid's coordinates of the unknowns
    std::vector<std::tuple<unsigned int, unsigned int>> v;

  public:
    /* Method to setup Ax=b system */
    int Setup(std::string inputfile);

    /* Method to create the solution vector structure v (map between grid point ij and its corresponding position in x solution vector). v = [(i,j), ....] ; gridpoint_to_x = {(1,0): 0; (2,0): 0, ....}. */
    std::map<std::tuple<int, int>, int> Create_GridVecMapping(std::vector<double> &x, std::vector<std::tuple<unsigned int, unsigned int>> &v, int &n_eqs, const int &ncols_input, const int &nrows_input);

    /* Method to fill the entries of matrix A in COO format and vector b according to the boundary conditions. */
    void Fill_Ab(std::vector<std::tuple<unsigned int, unsigned int>> &v, unsigned const int &n_eqs, const long unsigned int &ncols_input, std::vector<double> &b, const double &Th, const double &Tc, const double &length, const double &h);

    /* Method to compute the cool jet temperature for given column position on the grid.*/
    double Compute_Temp(const int grid_j, const double Tc, const double length, const double h);

    
    /* Method to solve system using CGsolver */
    int Solve(std::string soln_prefix);

};

#endif /* HEAT_HPP */