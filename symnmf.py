import sys
import copy
import math
import numpy as np
#import symnmfmodule

def algo(k, goal, vectors):
    # Get Dimensions
    d = len(vectors[0])
    n = len(vectors)

    # Call the relevant algorithm according to the given goal
    if goal == "symnmf":
        np.random.seed(1234)
        normalized_matrix = symnmfmodule.norm(vectors)
        total = sum(sum(row) for row in normalized_matrix)
        m = total/(n*n)
        initial_metrix = np.random.uniform(0, 2*math.sqrt(m/k), size=(n,k)) 
        result_matrix = symnmfmodule.symnmf(k, initial_metrix, normalized_matrix)
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
    for i in range(n):
        for j in range(k):
            if j==(k-1):
                print("{:.4f}".format(result_matrix[i][j]))
            else:
                print("{:.4f}".format(result_matrix[i][j])+",", end='')



def main():
    # Parse Arguments 
    if (len(sys.argv) == 4):
        k = sys.argv[1]
        goal = sys.argv[2]
        file_name = sys.argv[3]
    else:
        print("An Error Has Occurred")
        sys.exit(1)

    # Load files into NumPy arrays and convert it to lists of lists
    data = np.loadtxt(file_name, delimiter=',')
    X_datapoints = [datapoint.tolist() for datapoint in data]

    # Call the algo function of the relevant goal
    #algo(k, goal, X_datapoints)


main()