import sys
import copy
import math
import numpy as np
from sklearn.metrics import silhouette_score
import symnmfmodule

def kmeans_algo(k, iter, vectors):
    # Initialize
    epsilon = 0.001
    d = len(vectors[0])
    n = len(vectors)

    # Initialized centroids as first k vectors.
    centroids = []
    for i in range(k):
        centroids.append(copy.deepcopy(vectors[i]))

    # The k-means algorithm.
    for count_iter in range(iter):
        count_assign_vectors = [0 for _ in range(k)]
        sum_vectors_assign = [[0 for _ in range(d)] for _ in range(k)]
        labels = [0 for _ in range(n)]
        # Find for every vector the closest cetreoids.
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
            
            for index in range(d):
                sum_vectors_assign[min_center][index] += vectors[i][index]
            count_assign_vectors[min_center] += 1

        # Calculate the new cetreoids.
        prev_centroids = copy.deepcopy(centroids)
        for c in range(k):
            if (count_assign_vectors[c] != 0):
                for index in range(d):
                    centroids[c][index] = sum_vectors_assign[c][index]/count_assign_vectors[c]
        
        # Check if all centroids change less than epsilon in compare to there previous value.
        flag = True
        for i in range(k):
            sum = 0
            for j in range(d):
                delta_j = (centroids[i][j]-prev_centroids[i][j])**2
                sum += delta_j
            delta = math.sqrt(sum)
            if delta >= epsilon:
                flag = False
        if flag == True:
            break

    return labels


def main():
    # Parse Arguments 
    if (len(sys.argv) == 3):
        k = sys.argv[1]
        file_name = sys.argv[2]
    else:
        print("An Error Has Occurred")
        sys.exit(1)

    # Load files into NumPy arrays
    X_datapoints = np.loadtxt(file_name, delimiter=',')
    # Get Dimensions
    n, d = print(X_datapoints.shape)
    vectors = [datapoint.tolist() for datapoint in data]

    # Perform symNMF algorithm
    np.random.seed(1234)
    normalized_matrix = symnmfmodule.norm(vectors)
    total = sum(sum(row) for row in normalized_matrix)
    m = total/(n*n)
    H_initial_metrix = np.random.uniform(0, 2*math.sqrt(m/k), size=(n,n))
    H_final_matrix = symnmfmodule.symnmf(k, H_initial_metrix, normalized_matrix)

    # Find every row's cluster at H metrix.
    for i in range(k):
        labels_symnmf = [row.index(max(row)) for row in H_final_matrix]

    # Calculates symMNF score and Kmeans score.
    symnmf_score = silhouette_score(X_datapoints, labels_symnmf)
    print("nmf: {:.4f}".format(symnmf_score))
    labels_kmeans = kmeans_algo(3, 600, vectors)
    kmeans_score = silhouette_score(X_datapoints, labels_kmeans)
    print("kmeans: {:.4f}".format(kmeans_score))


main()