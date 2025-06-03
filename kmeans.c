#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef INFINITY
#define INFINITY (1.0/0.0)
#endif

#define EPSILON 0.001

struct cord {
    double value;
    struct cord *next;
};

struct vector {
    struct vector *next;
    struct cord *cords;
};

/* Checking whether the number is natural*/
int is_positive_integer(const char *str) {
    if (*str == '\0') return 0;  

    while (*str) {
        if (*str == '.') {
            str++;
            while (*str) {
                if (*str != '0') return 0;  
                str++;
            }
            return 1;

        } else if (*str < '0' || *str > '9') {
        return 0;
    }

        str++;
    }

    return 1;
}

/* Copying a linked list of coordinates */
struct cord* copy_cords(struct cord *source) {
    struct cord *head, *curr;
    if (!source) return NULL;

    head = malloc(sizeof(struct cord));
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

/* Print a vector */
void print_vector(struct cord *cords) {
    while (cords) {
        printf("%.4f", cords->value);
        if (cords->next) printf(",");
        cords = cords->next;
    }
    printf("\n");
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

int main(int argc, char **argv) {
    int k, iter = 400;
    int i, n = 0, count_iter, stop;
    int *counts;
    int min_index;
    double val, min_dist, delta;
    char ch;

    struct vector *head_vec, *curr_vec, *vec_iter;
    struct cord *head_cord, *curr_cord;
    struct cord **centroids;
    struct cord **sums;

    /* Input arguments check */
    if (argc != 2 && argc != 3) {
        printf("An Error Has Occurred\n");
        exit(1);
    }

    k = atoi(argv[1]);
    if (k <= 1) {
        printf("Incorrect number of clusters!\n");
        exit(1);
    }
       /* Checking whether k is a natural number*/
    if (!is_positive_integer(argv[1])) {
        printf("Incorrect number of clusters!\n");
        exit(1);
    }
    /* Checking whether iter is a natural number*/
    if (argc == 3) {
        iter = atoi(argv[2]);
        if (iter <= 1 || iter >= 1000) {
            printf("Incorrect maximum iteration!\n");
            exit(1);
        }
        if (!is_positive_integer(argv[2])) {
            printf("Incorrect maximum iteration!\n");
            exit(1);
        }
    }

    /* Reading vectors from input */
    head_vec = malloc(sizeof(struct vector));
    curr_vec = head_vec;
    curr_vec->next = NULL;

    head_cord = malloc(sizeof(struct cord));
    curr_cord = head_cord;
    curr_cord->next = NULL;

    while (scanf("%lf%c", &val, &ch) == 2) {
        curr_cord->value = val;
        if (ch == '\n') {
            curr_vec->cords = head_cord;
            curr_vec->next = malloc(sizeof(struct vector));
            curr_vec = curr_vec->next;
            curr_vec->next = NULL;
            curr_vec->cords = NULL;
            head_cord = malloc(sizeof(struct cord));
            curr_cord = head_cord;
            curr_cord->next = NULL;
            curr_cord->value = 0.0;
            n++;
            continue;
        } else {
            curr_cord->next = malloc(sizeof(struct cord));
            curr_cord = curr_cord->next;
            curr_cord->next = NULL;
            curr_cord->value = 0.0;
        }
    }       

    if (k >= n) {
        printf("Incorrect number of clusters!\n");
        free_vectors(head_vec);
        free(head_cord);
        return 1;
    }

    /* Initialize centroids with first k vectors */
    vec_iter = head_vec;
    centroids = calloc(k, sizeof(struct cord *));
    for (i = 0; i < k; i++) {
        centroids[i] = copy_cords(vec_iter->cords);
        vec_iter = vec_iter->next;
    }

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

    /* Output centroids */
    for (i = 0; i < k; i++) {
        print_vector(centroids[i]);
        free_cords(centroids[i]);
    }

    free(centroids);
    free_vectors(head_vec);
    free(head_cord);

    return 0;
}
