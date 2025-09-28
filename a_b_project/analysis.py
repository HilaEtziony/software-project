########################### IMPORTS ###########################

import sys
import copy
import math
import numpy as np
from sklearn.metrics import silhouette_score
import symnmfmodule
import shared

########################### CONSTS ###########################

EPSILON = 0.0001
MAX_ITER = 300

########################### FUNCTIONS ###########################
  
def is_converge(k, d, centroids, prev_centroids):
    """
    Check if all centroids change less than epsilon in compare to their previous value.
    :param k: number of clusters
    :type k: int
    :param d: the dimension of the vectors
    :type d: int
    :param centroids: current centroids
    :type centroids: list of lists
    :param prev_centroids: previous centroids
    :type prev_centroids: list of lists
    :return: True if all centroids change less than epsilon, False otherwise
    :rtype: bool
    """
    flag = True
    for i in range(k):
        sum = 0
        for j in range(d):
            delta_j = (centroids[i][j]-prev_centroids[i][j])**2
            sum += delta_j
        delta = math.sqrt(sum)
        if delta >= EPSILON:
            flag = False
    return flag

def calculate_new_centroids(k, d, centroids, sum_vectors_assign, count_assign_vectors):
    """
    Calculate the new centroids - for every cluster calculate the mean of the vectors assigned to it.
    :param k: number of clusters
    :type k: int
    :param d: the dimension of the vectors
    :type d: int
    :param centroids: current centroids
    :type centroids: list of lists
    :param sum_vectors_assign: the sum of the vectors assigned to each cluster for each element in the vector
    :type sum_vectors_assign: list of lists
    :param count_assign_vectors: the number of vectors assigned to each cluster
    :type count_assign_vectors: list
    :return: new centroids
    :rtype: list of lists
    """
    for cluster in range(k):
            if (count_assign_vectors[cluster] != 0):
                for index in range(d):
                    centroids[cluster][index] = sum_vectors_assign[cluster][index]/count_assign_vectors[cluster]
    return centroids

def init_centroids(k, vectors):
    """
    Initialize the centroids as the first k vectors.
    :param k: number of clusters
    :type k: int
    :param vectors: list of vectors
    :type vectors: list of lists
    :return: initialized centroids
    :rtype: list of lists
    """
    centroids = []
    for i in range(k):
        centroids.append(copy.deepcopy(vectors[i]))
    return centroids

def get_dimensions(vectors):
    """
    Get the dimensions of the vectors.
    :param vectors: list of vectors
    :type vectors: list of lists
    :return: n - number of vectors, d - dimension of the vectors
    :rtype: tuple (int, int)
    """
    n = len(vectors)
    d = len(vectors[0])
    return n, d

def kmeans_initialize_for_calc(k, d, n,):
    """
    Initialize the variables needed for the k-means calculation.
    :param k: number of clusters
    :type k: int
    :param d: the dimension of the vectors
    :type d: int
    :param n: number of vectors
    :type n: int
    :return: count_assign_vectors - the number of vectors assigned to each cluster,
             sum_vectors_assign - the sum of the vectors assigned to each cluster for each element in the vector,
             labels - the labels assigned to the vectors
    :rtype: tuple (list, list of lists, list)
    """
    count_assign_vectors = [0 for _ in range(k)]
    sum_vectors_assign = [[0 for _ in range(d)] for _ in range(k)]
    labels = [0 for _ in range(n)]
    return count_assign_vectors, sum_vectors_assign, labels

def kmeans_assign_labels(k, centroids, vectors):
    """
    The k-means algorithm - assign labels to the vectors according to the closest centroid.
    :param k: number of clusters
    :type k: int
    :param centroids: current centroids
    :type centroids: list of lists
    :param vectors: list of vectors
    :type vectors: list of lists
    :return: labels assigned to the vectors
    :rtype: list
    """
    n , d = get_dimensions(vectors)
    # The k-means algorithm.
    for count_iter in range(MAX_ITER):
        count_assign_vectors, sum_vectors_assign, labels = kmeans_initialize_for_calc(k, d, n)
        # Find for every vector the closest centroid.
        for i in range(n):
            min_delta = float('inf')
            min_center = 0
            for center_num in range(k):
                s = 0
                for j in range(d):
                    delta_j = (vectors[i][j]-centroids[center_num][j])**2
                    s += delta_j
                delta = math.sqrt(s)
                if delta < min_delta:
                    min_delta = delta
                    min_center = center_num
            labels[i] = min_center
            # Update the sum of the vectors assigned to the cluster and the count of the vectors assigned to it.
            for index in range(d):
                sum_vectors_assign[min_center][index] += vectors[i][index]
            count_assign_vectors[min_center] += 1
        # Calculate the new cetreoids.
        prev_centroids = copy.deepcopy(centroids)
        centroids = calculate_new_centroids(k, d, centroids, sum_vectors_assign, count_assign_vectors)
        # Check if the algorithm has converged.
        if is_converge(k, d, centroids, prev_centroids): break
    return labels

def kmeans_algorithm(k, vectors):
    """
    The k-means algorithm - assign labels to the vectors according to the closest centroid.
    :param k: number of clusters
    :type k: int
    :param vectors: list of vectors
    :type vectors: list of lists
    :return: labels assigned to the vectors
    :rtype: list
    """
    # Initialized centroids as first k vectors.
    centroids = init_centroids(k, vectors)
    # Assign labels to the vectors according to the closest centroid.
    labels = kmeans_assign_labels(k, centroids, vectors)
    return labels

def print_symNMF_score(k, X_datapoints):
    """
    Calculate and print the silhouette score of the symNMF algorithm.
    :param k: number of clusters
    :type k: int
    :param X_datapoints: data points
    :type X_datapoints: numpy.ndarray
    """
    # Perform symNMF algorithm
    normalized_matrix = symnmfmodule.norm(X_datapoints.tolist())
    H_initial_metrix = shared.initialize_H_Matrix(normalized_matrix, k)
    H_final_matrix = symnmfmodule.symnmf(H_initial_metrix, normalized_matrix)

    # Find every row's cluster at H metrix.
    labels_symnmf = [row.index(max(row)) for row in H_final_matrix]
    # Calculates symMNF score.
    symnmf_score = silhouette_score(X_datapoints, labels_symnmf)
    print("nmf: {:.4f}".format(symnmf_score))

def print_Kmeans_score(k, X_datapoints):
    """
    Calculate and print the silhouette score of the k-means algorithm.
    :param k: number of clusters
    :type k: int
    :param X_datapoints: data points
    :type X_datapoints: numpy.ndarray
    """
    # Perform k-means algorithm
    labels_kmeans = kmeans_algorithm(k, X_datapoints.tolist())
    # Calculates Kmeans score.
    kmeans_score = silhouette_score(X_datapoints, labels_kmeans)
    print("kmeans: {:.4f}".format(kmeans_score))

########################### MAIN ###########################

if __name__ == "__main__":
    try:
        # Parse Arguments 
        if (len(sys.argv) == 3):
            k = sys.argv[1]
            file_name = sys.argv[2]
        else:
            print("An Error Has Occurred")
            sys.exit(1)

        # Check if file's name is string and ends with .txt, exit if not valid
        shared.validate_filename(file_name)

        # Load file into NumPy arrays
        X_datapoints = np.loadtxt(file_name, delimiter=',', ndmin=2)

        # Get Dimensions
        n, d = X_datapoints.shape

        # Check if k is an integer and 1 < k < n, exit if not valid
        k = shared.check_and_get_k(k, n)
        
        # Print the silhouette scores of the symNMF and k-means algorithms.
        print_symNMF_score(k, X_datapoints)
        print_Kmeans_score(k, X_datapoints)

    except Exception:
        print("An Error Has Occurred")
        sys.exit(1)