import sys
import copy
import math
import pandas as pd
import numpy as np
import mykmeanspp

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

    # Reading files
    df1 = pd.read_csv(file_name_1, header=None)
    df2 = pd.read_csv(file_name_2, header=None)
    # inner join by key
    df_joined = pd.merge(df1, df2, on=0, how='inner')
    # sort by key
    df_joined = df_joined.sort_values(by=0)
    # convert to numpy array
    vectors = df_joined.to_numpy()
    # remove keys
    vectors = vectors[:, 1:]

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
    random_index = np.random.randint(1, len(vectors)) - 1
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
                    delta_j = (vectors[i,j]-centroids[center_num][j])**2
                    s += delta_j
                delta = math.sqrt(s)
                if delta < min_delta:
                    min_delta = delta
                    min_center = center_num
            prob.append(min_delta)
        # Find new cantroid.
        sum_prob = sum(prob)
        divided = [x / sum_prob for x in prob]
        sample = np.random.choice(np.arange(0, len(vectors)), size=1, p=divided)
        centroids.append(copy.deepcopy(vectors[sample[0]]))
        centroids_index.append(sample[0])

    centroids = mykmeanspp.fit(iter, epsilon, [arr.tolist() for arr in centroids], vectors.tolist())
  
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