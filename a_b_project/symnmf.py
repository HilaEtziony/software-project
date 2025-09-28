########################### IMPORTS ###########################

import sys
import numpy as np
import shared
import symnmfmodule

########################### FUNCTIONS ###########################

def print_matrix(matrix):
    """
    Print matrix in the required format.
    :param matrix: matrix for printing
    :type matrix: list of lists
    """
    for row in matrix:
        print(",".join(["{:.4f}".format(val) for val in row]))

def algorithm(k, goal, vectors):
    """
    Run the relevant algorithm according to the given goal.
    :param k: number of clusters
    :type k: int
    :param goal: the goal of the algorithm - "symnmf", "sym", "ddg", "norm"
    :type goal: str
    :param vectors: list of vectors
    :type vectors: list of lists
    """
    # Call the relevant algorithm according to the given goal
    if goal == "symnmf":
        normalized_matrix = symnmfmodule.norm(vectors)
        initial_metrix = shared.initialize_H_Matrix(normalized_matrix, k)
        result_matrix = symnmfmodule.symnmf(initial_metrix, normalized_matrix)
    elif goal == "sym":
        result_matrix = symnmfmodule.sym(vectors)
    elif goal == "ddg":
        result_matrix = symnmfmodule.ddg(vectors)
    elif goal == "norm":
        result_matrix = symnmfmodule.norm(vectors)
    else:
        print("An Error Has Occurred")
        sys.exit(1)

    # Print the metrix.
    print_matrix(result_matrix)

########################### MAIN ###########################

if __name__ == "__main__":
    try:
        # Parse Arguments 
        if (len(sys.argv) == 4):
            k = sys.argv[1]
            goal = sys.argv[2]
            file_name = sys.argv[3]
        else:
            print("An Error Has Occurred")
            sys.exit(1)
        
        # Check if file's name is string and ends with .txt, exit if not valid
        shared.validate_filename(file_name)

        # Load files into NumPy arrays
        data = np.loadtxt(file_name, delimiter=',', ndmin=2)

        # Get Dimensions
        n, d = data.shape

        # Check if k is an integer and 1 < k < n, exit if not valid
        k = shared.check_and_get_k(k, n)
        
        # Check if goal is valid, exit if not valid
        shared.validate_goal(goal)

        # Call the algorithm function with the relevant goal
        algorithm(k, goal, data.tolist())

    except Exception:
        print("An Error Has Occurred")
        sys.exit(1)