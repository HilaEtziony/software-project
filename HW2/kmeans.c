#include <stdio.h>
#include <stdlib.h>
#include "kmeans.h"
#include <math.h>

#ifndef INFINITY
#define INFINITY (1.0/0.0)
#endif


/* Copying a linked list of coordinates */
struct cord* copy_cords(struct cord *source) {
    struct cord *head, *curr;
    if (!source) return NULL;

    head = malloc(sizeof(struct cord));
    if (!head) return NULL;
    
    curr = head;
    curr->value = source->value;
    source = source->next;

    while (source) {
        curr->next = malloc(sizeof(struct cord));
        curr = curr->next;
        curr->value = source->value;
        source = source->next;
    }
    curr->next = NULL;
    return head;
}

/* Add one vector to another */
void add_vector(struct cord *sum, struct cord *v) {
    while (sum && v) {
        sum->value += v->value;
        sum = sum->next;
        v = v->next;
    }
}

/* Divide vector by a scalar */
void divide_vector(struct cord *v, int divisor) {
    while (v) {
        v->value /= divisor;
        v = v->next;
    }
}

/* Euclidean distance between vectors */
double delta_between_vectors(struct cord *a, struct cord *b) {
    double sum = 0.0, delta = 0.0;
    while (a && b) {
        delta = a->value - b->value;
        sum += delta * delta;
        a = a->next;
        b = b->next;
    }
    return sqrt(sum);
}


/* Free a list of cords */
void free_cords(struct cord *head) {
    struct cord *temp;
    while (head) {
        temp = head;
        head = head->next;
        free(temp);
    }
}

/* Free a list of vectors */
void free_vectors(struct vector *head) {
    struct vector *temp;
    while (head) {
        free_cords(head->cords);
        temp = head;
        head = head->next;
        free(temp);
    }
}
/* The main KMeans algorithm logic */
void kmeans_fit(int k, int iter, double epsilon, struct cord** centroids, struct vector* head_vec) {
    int i, count_iter, stop;
    int *counts;
    int min_index;
    double min_dist, delta;

    struct vector *vec_iter;
    struct cord **sums;

    /* Main K-Means loop */
    for (count_iter = 0; count_iter < iter; count_iter++) {
        counts = calloc(k, sizeof(int));
        sums = calloc(k, sizeof(struct cord *));
        for (i = 0; i < k; i++) {
            struct cord *curr, *ref;
            sums[i] = calloc(1, sizeof(struct cord));
            curr = sums[i];
            ref = centroids[i];
            while (ref->next) {
                curr->next = calloc(1, sizeof(struct cord));
                curr = curr->next;
                ref = ref->next;
            }
        }

        /* Assign vectors to nearest centroid */
        vec_iter = head_vec;
        while (vec_iter && vec_iter->cords) {
            min_dist = INFINITY;
            min_index = 0;

            for (i = 0; i < k; i++) {
                delta = delta_between_vectors(vec_iter->cords, centroids[i]);
                if (delta < min_dist) {
                    min_dist = delta;
                    min_index = i;
                }
            }

            add_vector(sums[min_index], vec_iter->cords);
            counts[min_index]++;
            vec_iter = vec_iter->next;
        }

        /* Recompute centroids */
        stop = 1;
        for (i = 0; i < k; i++) {
            if (counts[i] != 0)
                divide_vector(sums[i], counts[i]);

            delta = delta_between_vectors(centroids[i], sums[i]);
            if (delta >= epsilon)
                stop = 0;

            free_cords(centroids[i]);
            centroids[i] = copy_cords(sums[i]);
            free_cords(sums[i]);
        }

        free(sums);
        free(counts);

        if (stop)
            break;
    }
}
