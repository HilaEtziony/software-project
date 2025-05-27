#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define EPSILON 0.001

struct cord {
    double value;
    struct cord *next;
};

struct vector {
    struct vector *next;
    struct cord *cords;
};


// Copying a vector to another vector
struct cord* copy_cords(struct cord *source) {
    if (!source) return NULL;

    struct cord *head = malloc(sizeof(struct cord));
    struct cord *curr = head;
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


// Adding the vector value to sum
void add_vector(struct cord *sum, struct cord *v) {
    while (sum && v) {
        sum->value += v->value;
        sum = sum->next;
        v = v->next;
    }
}

// Calculating new centroid by dividing the sum of vectors by the number of vectors
void divide_vector(struct cord *v, int divisor) {
    while (v) {
        v->value /= divisor;
        v = v->next;
    }
}

// Calculating distance between vectors by Euclidean distance
double delta_between_vectors(struct cord *a, struct cord *b) {
    double sum = 0;
    while (a && b) {
        double delta = a->value - b->value;
        sum += delta * delta;
        a = a->next;
        b = b->next;
    }
    return sqrt(sum);
}

// Printing vector values
void print_vector(struct cord *cords) {
    while (cords) {
        printf("%.4f", cords->value);
        if (cords->next) printf(",");
        cords = cords->next;
    }
    printf("\n");
}

// Clearing memory starting from any struct head
void free_cords(struct cord *head) {
    struct cord *temp;
    while (head) {
        temp = head;
        head = head->next;
        free(temp);
    }
}
// Clearing memory starting from any vector head
void free_vectors(struct vector *head) {
    struct vector *temp;
    while (head) {
        free_cords(head->cords);
        temp = head;
        head = head->next;
        free(temp);
    }
}

int main(int argc, char **argv) {
    int k, iter = 400;
    int i, n = 0;

    if (argc != 2 && argc != 3) {
        printf("An Error Has Occurred\n");
        return 1;
    }

    k = atoi(argv[1]);
    if (k <= 1) {
        printf("Incorrect number of clusters!\n");
        return 1;
    }

    if (argc == 3) {
        iter = atoi(argv[2]);
        if (iter <= 1 || iter >= 1000) {
            printf("Incorrect maximum iteration!\n");
            return 1;
        }
    }

    // Process data from given file
    struct vector *head_vec = malloc(sizeof(struct vector));
    struct vector *curr_vec = head_vec;
    struct cord *head_cord = malloc(sizeof(struct cord));
    struct cord *curr_cord = head_cord;
    curr_cord->next = NULL;
    curr_vec->next = NULL; 

    double val;
    char ch;
    while (scanf("%lf%c", &val, &ch) == 2) {
        curr_cord->value = val;
        if (ch == '\n') {
            curr_vec->cords = head_cord;
            curr_vec->next = malloc(sizeof(struct vector));
            curr_vec = curr_vec->next;
            curr_vec->next = NULL;
            head_cord = malloc(sizeof(struct cord));
            curr_cord = head_cord;
            curr_cord->next = NULL;
            n++;
        } else {
            curr_cord->next = malloc(sizeof(struct cord));
            curr_cord = curr_cord->next;
            curr_cord->next = NULL;
        }
    }

    if (k >= n) {
        printf("Incorrect number of clusters!\n");
        return 1;
    }


    // Initialized centroids as first k vectors
    struct vector *vec_iter = head_vec;
    struct cord **centroids = malloc(k * sizeof(struct cord *));
    for (i = 0; i < k; i++) {
        centroids[i] = copy_cords(vec_iter->cords);
        vec_iter = vec_iter->next;
    }
    // The k-means algorithm.

    //Creating zeroing vectors for the scheme
    for (int count_iter = 0; count_iter < iter; count_iter++) {
        int *counts = calloc(k, sizeof(int));
        struct cord **sums = malloc(k * sizeof(struct cord *));
        for (i = 0; i < k; i++) {
            sums[i] = calloc(1, sizeof(struct cord));
            struct cord *curr = sums[i];
            struct cord *ref = centroids[i];
            while (ref->next) {
                curr->next = calloc(1, sizeof(struct cord));
                curr = curr->next;
                ref = ref->next;
            }
        }
        //Calculating the distance between any vector to his nearest centroid
        vec_iter = head_vec;
        while (vec_iter && vec_iter->cords) {
            double min_dist = INFINITY;
            int min_index = 0;
            for (i = 0; i < k; i++) {
                double dist = delta_between_vectors(vec_iter->cords, centroids[i]);
                if (dist < min_dist) {
                    min_dist = dist;
                    min_index = i;
                }
            }
            //update Sum and count after Calculating the nearest centroid
            add_vector(sums[min_index], vec_iter->cords);
            counts[min_index]++;
            vec_iter = vec_iter->next;
        }

        int stop = 1;
        for (i = 0; i < k; i++) {
            //New centroid calculation
            if (counts[i] != 0)
                divide_vector(sums[i], counts[i]);
            
            double delta = delta_between_vectors(centroids[i], sums[i]);
            if (delta >= EPSILON)
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

    for (i = 0; i < k; i++) {
        print_vector(centroids[i]);
        free_cords(centroids[i]);
    }

    free(centroids);
    free_vectors(head_vec);

    return 0;
}
