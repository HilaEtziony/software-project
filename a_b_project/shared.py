########################### IMPORTS ###########################

import numpy as np
import math

########################### FUNCTIONS ###########################

def initialize_H_Matrix(normalized_matrix, k):
    """
    Calculate initial H matrix for SymNMF algorithm using a random uniform distribution.
    :param normalized_matrix: normalized similarity matrix
    :type normalized_matrix: list of lists  
    :param k: number of clusters
    :type k: int
    :return: initial H matrix
    :rtype: list of lists
    """
    n = len(normalized_matrix)
    total = sum(sum(row) for row in normalized_matrix)
    m = total/(n*n)
    np.random.seed(1234)
    initial_metrix = np.random.uniform(0, 2*math.sqrt(m/k), size=(n,k)) 
    return initial_metrix.tolist()