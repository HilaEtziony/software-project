########################### IMPORTS ###########################

import numpy as np
import math
import sys

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

def string_is_integer(s):
    """
    Check if a string represents an integer.
    :param s: input string
    :type s: str
    :return: True if the string represents an integer, False otherwise
    :rtype: bool
    """
    try:
        return float(s).is_integer()
    except ValueError:
        return False

def check_and_get_k(k, n):
    """
    Check if k is an integer and 1 < k < n, and return k as an integer if valid, otherwise exit the program.
    :param n: a number to check k is smaller than it
    :type n: int  
    :param k: input k value as a string
    :type k: str    
    :return: k as an integer if valid
    :rtype: int
    """
    if (string_is_integer(k)):
        # Get k as an integer (also if it is given as a float string - for example 2.0)
        k = int(float(k))
        if not (1 < k and k < n):
            print("An Error Has Occurred")
            sys.exit(1)
        else:
            return k
    else:
        print("An Error Has Occurred")
        sys.exit(1)

def validate_filename(file_name):
    """
    Validate if the file name is a string and ends with .txt, otherwise exit the program
    :param file_name: input file name
    :type file_name: str
    """
    if not (isinstance(file_name, str) and file_name.endswith(".txt")):
        print("An Error Has Occurred")
        sys.exit(1)

def validate_goal(goal):
    """
    Validate if the goal is one of the accepted values, otherwise exit the program
    :param goal: input goal
    :type goal: str
    """
    if goal not in ["symnmf", "sym", "ddg", "norm"]:
        print("An Error Has Occurred")
        sys.exit(1)