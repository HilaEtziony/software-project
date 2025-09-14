#ifndef SYMNMF_H
#define SYMNMF_H

/* -------- data structures -------- */

struct cord {
    double       value;
    struct cord *next;
};

struct vector {
    struct cord   *cords;
    struct vector *next;
};

/* -------- exported helpers -------- */

void free_cords(struct cord *head);
void free_vectors(struct vector *head);
void free_dense(double **M, int n);
struct cord *cords_from_array(const double *arr, int m);

/* -------- exported algorithms -------- */

void symnmf_sym(struct cord **out_rows, struct vector *data_rows);
void symnmf_ddg(struct cord **out_rows, struct vector *data_rows);
void symnmf_norm(struct cord **out_rows, struct vector *data_rows);
void symnmf_symnmf(int k, struct cord **H, struct vector *W);

#endif /* SYMNMF_H */
