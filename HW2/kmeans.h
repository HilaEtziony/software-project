#ifndef KMEANS_H
#define KMEANS_H

struct cord {
    double value;
    struct cord *next;
};

struct vector {
    struct vector *next;
    struct cord *cords;
};

void kmeans_fit(int k, int iter, double epsilon, struct cord** centroids, struct vector* vectors);

void free_vectors(struct vector *head);

void free_cords(struct cord *head);
#endif
