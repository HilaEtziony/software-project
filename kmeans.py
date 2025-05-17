import sys
import copy
import math

def main():
    # Initialize
    epsilon = 0.001

    # Process data from given file to 2 dimensions's list
    lines = sys.stdin.readlines()
    vectors = [list(map(float, line.strip().split(','))) for line in lines]
    d = len(vectors[0])
    n = len(vectors)

    # Check k is valid
    arg1 = sys.argv[1]
    if (arg1.isdigit() and 1 < int(arg1) and int(arg1) < n):
        k = int(arg1)
    else:
        print("Incorrect number of clusters!")
        sys.exit(1)

    # Check iter is valid
    if (len(sys.argv) == 2): ## defualt iter
        iter = 400
    elif (len(sys.argv) == 3): 
        arg2 = sys.argv[2]
        if (arg2.isdigit() and 1 < int(arg2) and int(arg2) < 1000):
            iter = int(arg2)
        else:
            print("Incorrect maximum iteration!")
            sys.exit(1)
    else:
        print("An Error Has Occurred")
        sys.exit(1)

    # Initialized centroids as first k vectors.
    centroids = []
    for i in range(k):
        centroids.append(copy.deepcopy(vectors[i]))

    # The k-means algorithm.
    for count_iter in range(iter):
        count_assign_vectors = [0 for _ in range(k)]
        sum_vectors_assign = [[0 for _ in range(d)] for _ in range(k)]

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

    # Print the final k-means.
    for i in range(k):
        for j in range(d):
            if j==(d-1):
                print("{:.4f}".format(centroids[i][j]))
            else:
                print("{:.4f}".format(centroids[i][j])+",", end='')


main()