import sys
import copy
import math
import numpy as np
import mykmeanssp

def is_integer_string(s):
    try:
        return float(s).is_integer()
    except ValueError:
        return False
    
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def main():
    # Parse Arguments 
    if (len(sys.argv) == 5):
        k_before_validity = sys.argv[1]
        epsilon_before_validity = sys.argv[2]
        file_name_1 = sys.argv[3]
        file_name_2 = sys.argv[4]
    elif (len(sys.argv) == 6):
        k_before_validity = sys.argv[1]
        iter_before_validity = sys.argv[2]
        epsilon_before_validity = sys.argv[3]
        file_name_1 = sys.argv[4]
        file_name_2 = sys.argv[5]
    else:
        print("An Error Has Occurred")
        sys.exit(1)

    # Load files into NumPy arrays
    data1 = np.loadtxt(file_name_1, delimiter=',')
    data2 = np.loadtxt(file_name_2, delimiter=',')

    # Build dictionaries keyed by first column
    dict1 = {int(row[0]): row[1:] for row in data1}
    dict2 = {int(row[0]): row[1:] for row in data2}

    # Find common keys
    common_keys = sorted(set(dict1.keys()) & set(dict2.keys()))

    # Merge rows by key and build the list
    vectors = [np.concatenate([dict1[k], dict2[k]]).tolist() for k in common_keys]

    # Get Dimensions
    d = len(vectors[0])
    n = len(vectors)

    # Check k is valid
    if (is_integer_string(k_before_validity) and 1 < int(float(k_before_validity)) and int(float(k_before_validity)) < n):
        k = int(float(k_before_validity))
    else:
        print("Invalid number of clusters!")
        sys.exit(1)

    # Check iter is valid
    if (len(sys.argv) == 5): ## defualt iter
        iter = 300
    elif (len(sys.argv) == 6): 
        if (is_integer_string(iter_before_validity) and 1 < int(float(iter_before_validity)) and int(float(iter_before_validity)) < 1000):
            iter = int(float(iter_before_validity))
        else:
            print("Invalid maximum iteration!")
            sys.exit(1)
    else:
        print("An Error Has Occurred")
        sys.exit(1)

    # Check epsilon is valid
    if (is_number(epsilon_before_validity) and 0 <= float(epsilon_before_validity)):
        epsilon = float(epsilon_before_validity)
    else:
        print("Invalid epsilon!")
        sys.exit(1)

    # choose first centroid
    np.random.seed(1234)
    random_index = np.random.choice(len(vectors))
    centroids = [copy.deepcopy(vectors[random_index])]
    centroids_index = [random_index]

    # Find k centroids.
    for index in range(1,k):
        prob = []
        # Find for every vector the closest cetreoids.
        for i in range(len(vectors)):
            min_delta = float('inf')
            min_center = 0
            for center_num in range(len(centroids)):
                s = 0
                for j in range(d):
                    delta_j = (vectors[i][j]-centroids[center_num][j])**2
                    s += delta_j
                delta = math.sqrt(s)
                if delta < min_delta:
                    min_delta = delta
                    min_center = center_num
            prob.append(min_delta)
        # Find new cantroid.
        sum_prob = sum(prob)
        divided = [x / sum_prob for x in prob]
        divided = [float(x) for x in divided]
        sample = np.random.choice(len(vectors), p=divided)
        centroids.append(copy.deepcopy(vectors[sample]))
        centroids_index.append(sample)

    centroids = mykmeanssp.fit(iter, epsilon, centroids, vectors)

    # Print the centroids' original index.
    for j in range(k):
        if j==(k-1):
            print("{:}".format(centroids_index[j]))
        else:
            print("{:}".format(centroids_index[j])+",", end='')
    
    # Print the final k-means.
    for i in range(k):
        for j in range(d):
            if j==(d-1):
                print("{:.4f}".format(centroids[i][j]))
            else:
                print("{:.4f}".format(centroids[i][j])+",", end='')


main()