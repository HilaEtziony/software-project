from sklearn.datasets import load_iris
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt

def main():
    # Get data
    data = load_iris()

    k_values = list(range(1,11))
    inertia_values = []
    # Calculate inertia for each k
    for k in range(1,11):
        # Initialize k centroids
        centroids = KMeans(n_clusters=k, init='k-means++', random_state=0)
        # Run kmeans on data
        centroids.fit(data.data)
        # Calculate inetria
        inertia_values.append(centroids.inertia_)

    # Create points of k and its inetria
    points = []
    for i in range(10):
        points.append([k_values[i], inertia_values[i]])

    # Calculate the slope between the first and the last points
    main_slope = (points[9][1] - points[0][1]) / (points[9][0] - points[0][0])

    # Calculate the slope between every 2 adjacent point
    slopes = []
    for i in range(9):
        slopes.append((points[i+1][1] - points[i][1]) / (points[i+1][0] - points[i][0]))

    # Find the elbo
    for i in range(9):
        if slopes[i] <= main_slope and slopes[i+1] >= main_slope:
            elbow_index = i+1
            break
    elbow_k = k_values[elbow_index]

    # Create the graph
    plt.plot(k_values, inertia_values, marker='o')
    plt.xlabel('K ------>')
    plt.ylabel('Average Dispersion ------>')
    plt.title('Elbow Method for selection of optimal "K" clusters')
    plt.annotate(f'Elbow Point (k={elbow_k})',
             xy=(elbow_k, inertia_values[elbow_index]),
             xytext=(elbow_k + 0.5, inertia_values[elbow_index] + 20),
             arrowprops=dict(facecolor='black', arrowstyle='->'))
    plt.savefig('elbow.png')

main()