/* symnmf.c — algorithms + CLI (ANSI C / C90)
   Build: gcc -ansi -Wall -Wextra -Werror -pedantic-errors -O2 symnmf.c -o symnmf -lm
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "symnmf.h"

/* ---------------- tunables ---------------- */

#define EPS 1e-9
#define MAX_ITER 300
#define TOL 1e-4

/* ---------------- internal declarations ---------------- */

static int      cords_len(struct cord *c);
static int      vectors_len(struct vector *v);
static double** vectors_to_dense(struct vector *V, int *out_n, int *out_m);
static void     dense_to_cords_rows(double **M, int n, int m, struct cord **out_rows);

static void     die(void);
static char*    my_strdup(const char *s);
static struct vector* dense_to_vectors(double **M, int n, int m);
static struct vector* rows_to_vector(struct cord **rows, int n);
static void     print_rows(struct cord **rows, int n, int m);
static int      read_csv_dense(const char *path, double ***out_M, int *out_n, int *out_m);
static int      read_stdin_as_vectors(struct vector **out_head, int *out_n, int *out_m);

/* ---------------- memory helpers (exported in .h) ---------------- */

void free_cords(struct cord *head) {
    struct cord *tmp;
    tmp = NULL;
    while (head != NULL) {
        tmp = head;
        head = head->next;
        free(tmp);
    }
}

void free_vectors(struct vector *head) {
    struct vector *tmp;
    tmp = NULL;
    while (head != NULL) {
        tmp = head;
        head = head->next;
        free_cords(tmp->cords);
        free(tmp);
    }
}

/* ---------------- glue: linked lists <-> dense ---------------- */

static int cords_len(struct cord *c) {
    int len;
    len = 0;
    while (c != NULL) {
        len += 1;
        c = c->next;
    }
    return len;
}

static int vectors_len(struct vector *v) {
    int len;
    len = 0;
    while (v != NULL) {
        len += 1;
        v = v->next;
    }
    return len;
}

static double **vectors_to_dense(struct vector *V, int *out_n, int *out_m) {
    int i, j, n, m;
    struct vector *row;
    struct cord *c;
    double **M;

    i = 0; j = 0; n = 0; m = 0;
    row = NULL; c = NULL; M = NULL;

    if (V == NULL) {
        *out_n = 0;
        *out_m = 0;
        return NULL;
    }

    n = vectors_len(V);
    m = cords_len(V->cords);

    M = (double **)malloc(n * sizeof(*M));
    if (M == NULL) {
        *out_n = 0;
        *out_m = 0;
        return NULL;
    }

    for (i = 0; i < n; i += 1) {
        M[i] = (double *)malloc(m * sizeof(**M));
        if (M[i] == NULL) {
            for (j = 0; j < i; j += 1) {
                free(M[j]);
            }
            free(M);
            *out_n = 0;
            *out_m = 0;
            return NULL;
        }
    }

    row = V;
    for (i = 0; i < n && row != NULL; i += 1) {
        c = row->cords;
        for (j = 0; j < m && c != NULL; j += 1) {
            M[i][j] = c->value;
            c = c->next;
        }
        row = row->next;
    }

    *out_n = n;
    *out_m = m;
    return M;
}

/* exported in header */
struct cord *cords_from_array(const double *arr, int m) {
    int j;
    struct cord *head, *curr, *nnode;

    j = 0;
    head = NULL;
    curr = NULL;
    nnode = NULL;

    for (j = 0; j < m; j += 1) {
        nnode = (struct cord *)malloc(sizeof(struct cord));
        if (nnode == NULL) {
            free_cords(head);
            return NULL;
        }
        nnode->value = arr[j];
        nnode->next  = NULL;

        if (head == NULL) {
            head = nnode;
            curr = head;
        } else {
            curr->next = nnode;
            curr = curr->next;
            curr->next = NULL;
        }
    }
    return head;
}

static void dense_to_cords_rows(double **M, int n, int m, struct cord **out_rows) {
    int i;
    for (i = 0; i < n; i += 1) {
        out_rows[i] = cords_from_array(M[i], m);
    }
}

/* exported in header */
void free_dense(double **M, int n) {
    int i;
    if (M == NULL) {
        return;
    }
    for (i = 0; i < n; i += 1) {
        free(M[i]);
    }
    free(M);
}

/* ---------------- graph matrices: A, D, W (exported) ---------------- */

void sym(struct cord **out_rows, struct vector *data_rows) {
    int i, j, t, n, m, r;
    double **Xmat, **A;
    double s, df, v;

    i = j = t = n = m = r = 0;
    Xmat = NULL; A = NULL;
    s = 0.0; df = 0.0; v = 0.0;

    Xmat = vectors_to_dense(data_rows, &n, &m);
    if (n == 0 || m == 0) {
        return;
    }

    A = (double **)malloc(n * sizeof(*A));
    if (A == NULL) {
        free_dense(Xmat, n);
        return;
    }

    for (i = 0; i < n; i += 1) {
        A[i] = (double *)calloc(n, sizeof(**A));
        if (A[i] == NULL) {
            for (r = 0; r < i; r += 1) {
                free(A[r]);
            }
            free(A);
            free_dense(Xmat, n);
            return;
        }
    }

    /* A = exp(-||xi-xj||^2 / 2) */
    for (i = 0; i < n; i += 1) {
        for (j = i + 1; j < n; j += 1) {
            s = 0.0;
            for (t = 0; t < m; t += 1) {
                df = Xmat[i][t] - Xmat[j][t];
                s += df * df;
            }
            v = exp(-s / 2.0);
            A[i][j] = v;
            A[j][i] = v;
        }
    }

    dense_to_cords_rows(A, n, n, out_rows);

    free_dense(Xmat, n);
    free_dense(A, n);
}

/* Build D from X inside the function */
void ddg(struct cord **out_rows, struct vector *data_rows) {
    int i, j, n, m, r, n2;
    struct cord **A;
    struct vector *Avec;
    double **Amat, **D;
    double s;

    i = j = n = m = r = n2 = 0;
    A = NULL; Avec = NULL; Amat = NULL; D = NULL; s = 0.0;

    /* Build A(X) internally */
    n = vectors_len(data_rows);
    A = (struct cord **)calloc((size_t)n, sizeof(*A));
    if (A == NULL) {
        return;
    }
    sym(A, data_rows);

    /* Convert A to dense for easy row-sums */
    Avec = rows_to_vector(A, n);
    Amat = vectors_to_dense(Avec, &n2, &m);
    if (n2 == 0 || n2 != m) {
        free_vectors(Avec);
        free(A);
        free_dense(Amat, n2);
        return;
    }

    D = (double **)malloc(n * sizeof(*D));
    if (D == NULL) {
        free_dense(Amat, n);
        free_vectors(Avec);
        free(A);
        return;
    }
    for (i = 0; i < n; i += 1) {
        D[i] = (double *)calloc(n, sizeof(**D));
        if (D[i] == NULL) {
            for (r = 0; r < i; r += 1) {
                free(D[r]);
            }
            free(D);
            free_dense(Amat, n);
            free_vectors(Avec);
            free(A);
            return;
        }
    }
     /* D_ii = sum_j A_ij */
    for (i = 0; i < n; i += 1) {
        s = 0.0;
        for (j = 0; j < n; j += 1) {
            s += Amat[i][j];
        }
        D[i][i] = s;
    }

    dense_to_cords_rows(D, n, n, out_rows);

    free_dense(Amat, n);
    free_vectors(Avec);   /* frees the cords that Arows points to */
    free(A);          /* free only the pointer array */
    for (i = 0; i < n; i += 1) {
        free(D[i]);
    }
    free(D);
}

/* Build W = D^{-1/2} A D^{-1/2} from X inside the function */
void norm(struct cord **out_rows, struct vector *data_rows) {
    int i, j, n, m, r, n2;
    struct cord **Arows;
    struct vector *Avec;
    double **Amat, **W, *deg;
    double di, dj, s;

    i = j = n = m = r = n2 = 0;
    Arows = NULL; Avec = NULL; Amat = NULL; W = NULL; deg = NULL;
    di = 0.0; dj = 0.0; s = 0.0;

    /* Build A(X) internally */
    n = vectors_len(data_rows);
    Arows = (struct cord **)calloc((size_t)n, sizeof(*Arows));
    if (Arows == NULL) {
        return;
    }
    sym(Arows, data_rows);

    /* Convert to dense */
    Avec = rows_to_vector(Arows, n);
    Amat = vectors_to_dense(Avec, &n2, &m);
    if (n2 == 0 || n2 != m) {
        free_vectors(Avec);
        free(Arows);
        free_dense(Amat, n2);
        return;
    }

    deg = (double *)calloc(n, sizeof(*deg));
    if (deg == NULL) {
        free_dense(Amat, n);
        free_vectors(Avec);
        free(Arows);
        return;
    }
    for (i = 0; i < n; i += 1) {
        s = 0.0;
        for (j = 0; j < n; j += 1) {
            s += Amat[i][j];
        }
        deg[i] = s;
    }

    W = (double **)malloc(n * sizeof(*W));
    if (W == NULL) {
        free(deg);
        free_dense(Amat, n);
        free_vectors(Avec);
        free(Arows);
        return;
    }
    for (i = 0; i < n; i += 1) {
        W[i] = (double *)malloc(n * sizeof(**W));
        if (W[i] == NULL) {
            for (r = 0; r < i; r += 1) {
                free(W[r]);
            }
            free(W);
            free(deg);
            free_dense(Amat, n);
            free_vectors(Avec);
            free(Arows);
            return;
        }
    }
    /* W = D^{-1/2} A D^{-1/2} */
    for (i = 0; i < n; i += 1) {
        di = (deg[i] > 0.0) ? 1.0 / sqrt(deg[i]) : 0.0;
        for (j = 0; j < n; j += 1) {
            dj = (deg[j] > 0.0) ? 1.0 / sqrt(deg[j]) : 0.0;
            W[i][j] = di * Amat[i][j] * dj;
        }
    }

    dense_to_cords_rows(W, n, n, out_rows);

    free(deg);
    free_dense(Amat, n);
    free_vectors(Avec);
    free(Arows);
    for (i = 0; i < n; i += 1) {
        free(W[i]);
    }
    free(W);
}

/* ---------------- SymNMF multiplicative updates (exported) ---------------- */

void symnmf_algo(int k, struct cord **H, struct vector *W) {
    /* Multiplicative SymNMF with beta=0.5:
       H <- H * [ (1 - beta) + beta * ( (W H) ./ ((H H^T) H) ) ]
       Stops when L1 change across H is below TOL or after MAX_ITER.
       Assumes: W is n x n, H is array of n linked-lists each of length k. */
    int i, j, r, ii, jj, iter, rr, n_w, m_w, n;
    double **Wmat, **Hmat, **WH, **HHt, **HHTH;
    double Wij, s, denom, oldv, newv, diff_sum, beta;
    struct cord *node, *tail;

    i = j = r = ii = jj = iter = rr = n_w = m_w = n = 0;
    Wmat = Hmat = WH = HHt = HHTH = NULL;
    Wij = s = denom = oldv = newv = diff_sum = 0.0;
    node = tail = NULL;
    beta = 0.5; /* per specification */

    /* W: vector-of-rows (linked lists) -> dense and validate square */
    Wmat = vectors_to_dense(W, &n_w, &m_w);
    if (n_w == 0 || n_w != m_w) {
        free_dense(Wmat, n_w);
        return;
    }
    n = n_w;

    /* Build dense H (n x k) from the linked-lists in H */
    Hmat = (double **)malloc(n * sizeof(*Hmat));
    if (Hmat == NULL) { free_dense(Wmat, n); return; }
    for (i = 0; i < n; i += 1) {
        Hmat[i] = (double *)calloc((size_t)k, sizeof(**Hmat));
        if (Hmat[i] == NULL) {
            for (rr = 0; rr < i; rr += 1) free(Hmat[rr]);
            free(Hmat); free_dense(Wmat, n); return;
        }
        node = H[i];
        for (j = 0; j < k && node != NULL; j += 1) {
            Hmat[i][j] = (node->value >= 0.0) ? node->value : 0.0; /* keep non-negative */
            node = node->next;
        }
    }

       /* Allocate intermediate matrices:
       - WH    = W * H        (size n × k)   -> numerator in the update
       - HHt   = H * H^T      (size n × n)   -> Gram matrix of H
       - HHTH  = (H H^T) * H  (size n × k)   -> denominator in the update
    */
    WH   = (double **)malloc(n * sizeof(*WH));
    HHt  = (double **)malloc(n * sizeof(*HHt));
    HHTH = (double **)malloc(n * sizeof(*HHTH));
    if (WH == NULL || HHt == NULL || HHTH == NULL) {
        free(WH); free(HHt); free(HHTH);
        free_dense(Hmat, n); free_dense(Wmat, n);
        return;
    }
    for (i = 0; i < n; i += 1) {
        WH[i]   = (double *)calloc((size_t)k, sizeof(**WH));
        HHt[i]  = (double *)calloc((size_t)n, sizeof(**HHt));
        HHTH[i] = (double *)calloc((size_t)k, sizeof(**HHTH));
        if (WH[i] == NULL || HHt[i] == NULL || HHTH[i] == NULL) {
            int t;
            for (t = 0; t <= i; t += 1) { free(WH[t]); free(HHt[t]); free(HHTH[t]); }
            free(WH); free(HHt); free(HHTH);
            free_dense(Hmat, n); free_dense(Wmat, n);
            return;
        }
    }

    /* Main loop */
    for (iter = 0; iter < MAX_ITER; iter += 1) {

    /* ---------- WH = W * H   (n × k) ---------- */
    for (i = 0; i < n; i += 1) {
        for (r = 0; r < k; r += 1) {
            WH[i][r] = 0.0;
        }
        for (j = 0; j < n; j += 1) {
            Wij = Wmat[i][j];
            if (Wij != 0.0) {
                for (r = 0; r < k; r += 1) {
                    WH[i][r] += Wij * Hmat[j][r];
                }
            }
        }
    }

    /* ---------- HHt = H * H^T   (n × n) ---------- */
    for (ii = 0; ii < n; ii += 1) {
        for (jj = 0; jj < n; jj += 1) {
            s = 0.0;
            for (r = 0; r < k; r += 1) {
                s += Hmat[ii][r] * Hmat[jj][r];
            }
            HHt[ii][jj] = s;
        }
    }

    /* ---------- HHTH = (H H^T) * H   (n × k) ---------- */
    for (i = 0; i < n; i += 1) {
        for (r = 0; r < k; r += 1) {
            s = 0.0;
            for (j = 0; j < n; j += 1) {
                s += HHt[i][j] * Hmat[j][r];
            }
            HHTH[i][r] = s;
        }
    }

    /* ---------- Update step ---------- 
       H = H * ( (1 - beta) + beta * (WH ./ max(HHTH, EPS)) )
    */
    diff_sum = 0.0;
    for (i = 0; i < n; i += 1) {
        for (r = 0; r < k; r += 1) {
            denom = HHTH[i][r];
            if (denom < EPS) {
                denom = EPS;
            }

            oldv = Hmat[i][r];
            newv = oldv * ( (1.0 - beta) + beta * (WH[i][r] / denom) );

            if (newv < 0.0) {
                newv = 0.0; /* numerical safety */
            }

            Hmat[i][r] = newv;
            diff_sum += (newv - oldv)*(newv - oldv);
        }
    }

    if (diff_sum < TOL) {
        break;
    }
}


    /* Write back to the linked lists H (extend rows if needed) */
    for (i = 0; i < n; i += 1) {
        node = H[i]; tail = NULL;
        if (node != NULL) { while (node->next != NULL) node = node->next; tail = node; }
        node = H[i];
        for (r = 0; r < k; r += 1) {
            if (node != NULL) {
                node->value = Hmat[i][r];
                node = node->next;
            } else {
                struct cord *nc = (struct cord *)malloc(sizeof(struct cord));
                if (nc == NULL) break; /* best-effort */
                nc->value = Hmat[i][r];
                nc->next  = NULL;
                if (H[i] == NULL) { H[i] = nc; tail = nc; }
                else { tail->next = nc; tail = nc; }
            }
        }
    }

    /* Free */
    free_dense(WH, n);
    free_dense(HHt, n);
    free_dense(HHTH, n);
    free_dense(Hmat, n);
    free_dense(Wmat, n);
}




/* ---------------- CLI utilities ---------------- */

static void die(void) {
    /* Required exact error message for the automatic tests */
    fprintf(stderr, "An Error Has Occurred\n");
    exit(1);
}

static char *my_strdup(const char *s) {
    size_t n;
    char *p;
    n = strlen(s) + 1;
    p = (char *)malloc(n);
    if (p == NULL) {
        return NULL;
    }
    memcpy(p, s, n);
    return p;
}

static struct vector *dense_to_vectors(double **M, int n, int m) {
    int i;
    struct vector *head, *curr, *nv;

    head = NULL; curr = NULL; nv = NULL;

    for (i = 0; i < n; i += 1) {
        nv = (struct vector *)malloc(sizeof(struct vector));
        if (nv == NULL) {
            die();
        }
        nv->cords = cords_from_array(M[i], m);
        nv->next  = NULL;

        if (head == NULL) {
            head = nv;
            curr = nv;
        } else {
            curr->next = nv;
            curr = nv;
        }
    }
    return head;
}

static struct vector *rows_to_vector(struct cord **rows, int n) {
    int i;
    struct vector *head, *curr, *nv;

    head = NULL; curr = NULL; nv = NULL;

    for (i = 0; i < n; i += 1) {
        nv = (struct vector *)malloc(sizeof(struct vector));
        if (nv == NULL) {
            die();
        }
        nv->cords = rows[i];
        nv->next  = NULL;

        if (head == NULL) {
            head = nv;
            curr = nv;
        } else {
            curr->next = nv;
            curr = nv;
            curr->next = NULL;
        }
    }
    return head;
}

static void print_rows(struct cord **rows, int n, int m) {
    int i, j;
    struct cord *c;

    i = j = 0;
    c = NULL;

    for (i = 0; i < n; i += 1) {
        c = rows[i];
        for (j = 0; j < m; j += 1) {
            if (c == NULL) {
                die();
            }
            if (j == m - 1) {
                printf("%.4f\n", c->value);
            } else {
                printf("%.4f,", c->value);
            }
            c = c->next;
        }
    }
}

static int read_csv_dense(const char *path, double ***out_M, int *out_n, int *out_m) {
    FILE *f;
    const size_t BUFSZ = 1u << 20; /* 1MB */
    char *buf, *p, *tok, *line;
    int cap, n, m, j;
    int only_ws, cnt;
    int i;
    double **M;
    double **tmp;

    f = fopen(path, "r");
    if (f == NULL) {
        return -1;
    }

    cap = 128; n = 0; m = -1; j = 0;
    buf = (char *)malloc(BUFSZ);
    if (buf == NULL) {
        fclose(f);
        return -1;
    }

    M = (double **)malloc(cap * sizeof(*M));
    if (M == NULL) {
        free(buf);
        fclose(f);
        return -1;
    }

    while (fgets(buf, (int)BUFSZ, f) != NULL) {
        only_ws = 1;
        p = buf;
        while (*p != '\0') {
            if (!(*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n')) {
                only_ws = 0;
                break;
            }
            p += 1;
        }
        if (only_ws == 1) {
            continue;
        }

        if (m < 0) {
            cnt = 1;
            for (p = buf; *p != '\0'; p += 1) {
                if (*p == ',') {
                    cnt += 1;
                }
            }
            m = cnt;
        }

        line = my_strdup(buf);
        if (line == NULL) {
            free(buf);
            fclose(f);
            return -1;
        }

        j = 0;
        M[n] = (double *)malloc(m * sizeof(**M));
        if (M[n] == NULL) {
            free(line);
            free(buf);
            fclose(f);
            return -1;
        }

        tok = strtok(line, ",\r\n");
        while (tok != NULL && j < m) {
            M[n][j] = strtod(tok, NULL);
            j += 1;
            tok = strtok(NULL, ",\r\n");
        }
        free(line);

        if (j != m) {
            free(M[n]);
            free(buf);
            fclose(f);
            return -1;
        }

        n += 1;
        if (n == cap) {
            cap *= 2;
            tmp = (double **)realloc(M, cap * sizeof(*M));
            if (tmp == NULL) {
                free(buf);
                fclose(f);
                return -1;
            }
            M = tmp;
        }
    }

    free(buf);
    fclose(f);

    if (m <= 0 || n <= 0) {
        for (i = 0; i < n; i += 1) {
            free(M[i]);
        }
        free(M);
        return -1;
    }

    *out_M = M;
    *out_n = n;
    *out_m = m;
    return 0;
}

/* read from stdin: comma-separated values, newline ends row */
static int read_stdin_as_vectors(struct vector **out_head, int *out_n, int *out_m) {
    struct vector *head_vec, *curr_vec, *prev_vec;
    struct cord   *head_cord, *curr_cord;
    int n;
    double val;
    char ch;

    head_vec = (struct vector *)malloc(sizeof(struct vector));
    if (!head_vec) return -1;
    curr_vec = head_vec;
    curr_vec->next  = NULL;
    curr_vec->cords = NULL;
    prev_vec = NULL;

    head_cord = (struct cord *)malloc(sizeof(struct cord));
    if (!head_cord) { free(head_vec); return -1; }
    curr_cord = head_cord;
    curr_cord->next  = NULL;
    curr_cord->value = 0.0;

    n = 0;

    while (scanf("%lf%c", &val, &ch) == 2) {
        curr_cord->value = val;

        if (ch == '\n') {
            curr_vec->cords = head_cord;
            prev_vec = curr_vec;

            curr_vec->next = (struct vector *)malloc(sizeof(struct vector));
            if (!curr_vec->next) { free_vectors(head_vec); return -1; }
            curr_vec = curr_vec->next;
            curr_vec->next  = NULL;
            curr_vec->cords = NULL;

            head_cord = (struct cord *)malloc(sizeof(struct cord));
            if (!head_cord) { free_vectors(head_vec); return -1; }
            curr_cord = head_cord;
            curr_cord->next  = NULL;
            curr_cord->value = 0.0;

            n++;
        } else {
            curr_cord->next = (struct cord *)malloc(sizeof(struct cord));
            if (!curr_cord->next) { free_vectors(head_vec); return -1; }
            curr_cord = curr_cord->next;
            curr_cord->next  = NULL;
            curr_cord->value = 0.0;
        }
    }

    if (n == 0) {
        free(head_cord);
        free(head_vec);
        return -1;
    }

    free(head_cord);
    free(curr_vec);
    prev_vec->next = NULL;

    *out_head = head_vec;
    *out_n = n;
    *out_m = cords_len(head_vec->cords);
    return 0;
}

/* ---------------- main (CLI) ---------------- */

int main(int argc, char **argv) {
    int n, d, from_stdin, i;
    const char *goal, *path;
    struct vector *X_vecs, *Avec;
    double **X_dense;
    struct cord **A, **D, **Wrows;

    n = 0; d = 0; from_stdin = 0; i = 0;
    goal = NULL; path = NULL;
    X_vecs = NULL; X_dense = NULL;
    Avec = NULL; A = NULL; D = NULL; Wrows = NULL;

    if (argc != 3) {
        die();
    }

    goal = argv[1];
    path = argv[2];

    /* SymNMF is not supported in the CLI because k is not provided */
    if (strcmp(goal, "symnmf") == 0) {
        die();
    }

    /* Load input: CSV file or stdin ("-") */
    if (strcmp(path, "-") == 0) {
        from_stdin = 1;
        if (read_stdin_as_vectors(&X_vecs, &n, &d) != 0) {
            die();
        }
    } else {
        if (read_csv_dense(path, &X_dense, &n, &d) != 0) {
            die();
        }
        X_vecs = dense_to_vectors(X_dense, n, d);
    }

    if (strcmp(goal, "sym") == 0) {
        A = (struct cord **)calloc((size_t)n, sizeof(*A));
        if (A == NULL) {
            die();
        }
        sym(A, X_vecs);
        print_rows(A, n, n);
        for (i = 0; i < n; i += 1) {
            free_cords(A[i]);
        }
        free(A);

    } else if (strcmp(goal, "ddg") == 0) {
        D = (struct cord **)calloc((size_t)n, sizeof(*D));
        if (D == NULL) {
            die();
        }
        ddg(D, X_vecs);               /* pass X directly */
        print_rows(D, n, n);
        for (i = 0; i < n; i += 1) {
            free_cords(D[i]);
        }
        free(D);

    } else if (strcmp(goal, "norm") == 0) {
        Wrows = (struct cord **)calloc((size_t)n, sizeof(*Wrows));
        if (Wrows == NULL) {
            die();
        }
        norm(Wrows, X_vecs);          /* pass X directly */
        print_rows(Wrows, n, n);
        for (i = 0; i < n; i += 1) {
            free_cords(Wrows[i]);
        }
        free(Wrows);

    } else {
        die();
    }

    free_vectors(X_vecs);
    if (!from_stdin && X_dense != NULL) {
        free_dense(X_dense, n);
    }
    return 0;
}
