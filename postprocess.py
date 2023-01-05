import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse
import scipy.sparse.linalg
import warnings
import matplotlib
import sys
import math
import bisect
from scipy.interpolate import interp1d


def Get_GeometricParameters(input_file_list):
    """
        Receives a list of each line in the input file and get the parameters that define the geometry/system to analyze. Returns: length, width, h (spacing between points), Tc (cool jet temperature), and Th (hot temperature).
    """
    geometry_param = input_file_list[0].split() # Geometry parameters
    length = float(geometry_param[0])
    width = float(geometry_param[1])
    h = float(geometry_param[2])
    isoT_data = input_file_list[1].split() # Isothermal boundaries data
    Tc = int(isoT_data[0]) # Cool jet temperature
    Th = int(isoT_data[1]) # Hot boundary temperature

    return length, width, h, Tc, Th

def Create_Grid(nrows, ncols):
    """
        Receives the number of nrows and columns needed to create a grid/array for the system. Fills each entry with the corresponding solution's temperature data and the known temperature at the boundaries (Tc and Th from the input file).
    """
    grid = np.zeros((nrows, ncols)) # Create np array of nrows x ncols
    grid[0,:] = Th # Set temperature for lower boundary (cool jet)


    plt.plot(grid[0,:])
    plt.savefig("Last row of the distribution.png")

    grid[nrows-1,:] = Tc # Set temperature for upper boundary

    pos_val = 0 # Position of the value in the solution's vector
    for j in range(ncols-1):
        for i in range(nrows-2):
            grid[i+1, j] = float(sol_components[pos_val]) # Fill grid entry with temperature value
            pos_val += 1 # Next position
    grid[1:nrows-1, ncols-1] = grid[1:nrows-1, 0] # Fill the right boundary of the grid

    return grid

def Plot_Solution(grid, length, width, ncols, nrows):
    """
        Computes the temperature isoline for each column on the grid/array. Receives the grid/array, number of columns and the mean temperature. 

        Create a colormesh plot of the temperature distribution for the given geometry and solution file. Receives: length, width, ncols, and nrows. Saves the file with the name "Temperature_Distribution.png" on the working directory.
    """
    
    # T_isoline = []
    # for j in range(ncols):
    #     T_mean_col = np.mean(grid[:,j])
    #     T_isoline.append(T_mean_col)
    # x_pos = []
    # y_val = []
    # for i in range(len(T_isoline)):
    #         x_pos.append(i)
    #         y_val.append(T_isoline[i])
    #         T = T_isoline[i]
    # f = interp1d(x_pos,y_val,kind='linear')
    # xnew = np.arange(0,length,h)
    # ynew= f(xnew)
    # plt.plot(xnew,ynew,color = "black",lw=3)
    # for i in range(0,len(T_isoline)):
    #     # print(T_isoline[i])
    #     # print(T)
    #     if T_isoline[i] != T:
    #         x_pos.append(i)
    #         y_val.append(T_isoline[i])
    #         T = T_isoline[i]
    # x_pos.append(len(T_isoline))
    # y_val.insert(len(y_val)//2,y_val[len(y_val)//2])
    # x2 = np.arange(0,ncols,h)
    # f = interp1d(x_pos,y_val,kind='cubic')
    # y2= f(x2)


    X = np.linspace(0, length, ncols)
    Y = np.linspace(width, 0, nrows)
    # plt.plot(x2,y2,color = "black",lw=2.5)
    
    plt.pcolormesh(X,Y, grid)
    plt.ylim([(-1)*width, 2*width])
    plt.xlim([0, length])
    plt.colorbar()
    plt.savefig("Temperature_Distribution.png")




if __name__ == "__main__":
	if len(sys.argv) <= 2:
		# missing arguments , print usage message 
		print("Usage:")
		print(" $ python3 postprocess.py <input#.txt> <solution#.txt>")
		sys.exit(0)

input_file = open(sys.argv[1], "r") # Set input file on 'read' mode
solution_file = open(sys.argv[2], "r") # Set solutions file on 'read' mode
print("Input file processed: "+sys.argv[1])
# Read the content of the files
input_file_list = input_file.readlines() 
sol_components = solution_file.readlines() 
# Get geometric parameters defining the system
length, width, h, Tc, Th = Get_GeometricParameters(input_file_list)
# Number of nrows and columns in the matrix/grid
nrows = int(width/h)
ncols = int(length/h)
# Create and fill an array with the solution's temperature values
grid = Create_Grid(nrows, ncols)
# Compute the mean temperature
T_mean = np.mean(grid)
print("Mean Temperature: %.5f"%T_mean)
# Plot and save the solution
Plot_Solution(grid, length, width, ncols, nrows)
