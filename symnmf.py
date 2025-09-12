import sys
import copy
import math
import numpy as np
import symnmfmodule

def print_matrix(matrix):
    for row in matrix:
        print(",".join(["{:.4f}".format(val) for val in row]))

def algo(k, goal, data):
    # Get Dimensions
    n, d = data.shape
    vectors = data.tolist()

    # Call the relevant algorithm according to the given goal
    if goal == "symnmf":
        normalized_matrix = symnmfmodule.py_norm(vectors)
        # Calculate initial H
        total = sum(sum(row) for row in normalized_matrix)
        m = total/(n*n)
        np.random.seed(1234)
        initial_metrix = np.random.uniform(0, 2*math.sqrt(m/k), size=(n,k)) 
        initial_metrix =initial_metrix.tolist()
        result_matrix = symnmfmodule.py_symnmf(k, initial_metrix, normalized_matrix)
    elif goal == "sym":
        result_matrix = symnmfmodule.py_sym(vectors)
    elif goal == "ddg":
        result_matrix = symnmfmodule.py_ddg(vectors)
    elif goal == "norm":
        result_matrix = symnmfmodule.py_norm(vectors)
    else:
        print("An Error Has Occurred")
        sys.exit(1)

    # Print the metrix.
    print_matrix(result_matrix)

def main():
    # Parse Arguments 
    if (len(sys.argv) == 4):
        k = int(sys.argv[1])
        goal = sys.argv[2]
        file_name = sys.argv[3]
    else:
        print("An Error Has Occurred")
        sys.exit(1)

    # Load files into NumPy arrays
    data = np.loadtxt(file_name, delimiter=',', ndmin=2)

    # Call the algo function with the relevant goal
    algo(k, goal, data)


main()