#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include "symnmf.h"

#ifndef SIZE_MAX
#define SIZE_MAX ((size_t)-1)
#endif

#define TINY_POSITIVE 1e-9
#define MAX_ITER 300
#define EPS 1e-4


static void *xmalloc_size(size_t count, size_t elem) {
    if (count != 0 && elem > SIZE_MAX / count) {
        return NULL; /* overflow */
    }
    return malloc(count * elem);
}

static void *xcalloc_size(size_t count, size_t elem) {
    if (count != 0 && elem > SIZE_MAX / count) {
        return NULL; /* overflow */
    }
    return calloc(count, elem);
}

/* Count cords (elements) in a cord linked list. */
int cords_len(struct cord *c) {
    int len;
    len = 0;
    while (c != NULL) { 
        len += 1;
        c = c->next; 
     }
    return len;
}

/* Count vectors in a vector linked list. */
static int vectors_len(struct vector *v) {
    int len;
    len = 0;
    while (v != NULL) {
        len += 1;
        v = v->next;
    }
    return len;
}

/* Print error and exit */
static void die(void) {
    printf("An Error Has Occurred\n");
    exit(1);
}

/* ANSI-safe strdup: returns newly allocated copy or NULL. */
static char *my_strdup(const char *s) {
    size_t n;
    char *p;
    n = strlen(s) + 1;
    p = (char *)malloc(n);
    if (p == NULL) return NULL;
    memcpy(p, s, n);
    return p;
}

/* Build a cord list from a double row of length m. (exported in .h) */
struct cord *cords_from_array(const double *arr, int m) {
    int j;
    struct cord *head, *tail, *new_node;

    j = 0;
    head = NULL;
    tail = NULL;
    new_node = NULL;

    for (j = 0; j < m; j += 1) {
        new_node = (struct cord *)malloc(sizeof(struct cord));
        if (new_node == NULL) {
            free_cords(head);
            return NULL;
        }
        new_node->value = arr[j];
        new_node->next = NULL;
        if (head == NULL) {
            head = new_node;
            tail = head;
        } else {
            tail->next = new_node;
            tail = tail->next;
            tail->next = NULL;
        }
    }
    return head;
}

/* Convert dense M[n][m] to vector linked list. */
static struct vector *dense_to_vectors(double **M, int n, int m) {
    int i;
    struct vector *head, *tail, *new_vec;
    head = NULL; tail = NULL; new_vec = NULL;

    for (i = 0; i < n; i += 1) {
        new_vec = (struct vector *)malloc(sizeof(struct vector));
        if (new_vec == NULL) die();
        new_vec->cords = cords_from_array(M[i], m);
        new_vec->next = NULL;
        if (head == NULL) {
            head = new_vec;
            tail = new_vec;
        } else {
            tail->next = new_vec;
            tail = new_vec;
        }
    }
    return head;
}

/* Allocate dense M with n rows and m cols (zeroed). */
static double **alloc_dense(int n, int m) {
    int i, r;
    double **M;
    M = (double **)xmalloc_size((size_t)n, sizeof(*M));
    if (M == NULL) return NULL;
    for (i = 0; i < n; i += 1) {
        M[i] = (double *)xcalloc_size((size_t)m, sizeof(**M));
        if (M[i] == NULL) {
            for (r = 0; r < i; r += 1) free(M[r]);
            free(M);
            return NULL;
        }
    }
    return M;
}

/* Copy vector list V into a dense matrix, returning M and its n,m. */
static double **vectors_to_dense(struct vector *V, int *out_n, int *out_m) {
    int i, j, n, m;
    struct vector *row_vec;
    struct cord *cord_node;
    double **M;

    if (V == NULL) { *out_n = 0; *out_m = 0; return NULL; }
    n = vectors_len(V);
    m = cords_len(V->cords);
    M = alloc_dense(n, m);
    if (M == NULL) { *out_n = 0; *out_m = 0; return NULL; }
    row_vec = V;
    for (i = 0; i < n && row_vec != NULL; i += 1) {
        cord_node = row_vec->cords;
        for (j = 0; j < m && cord_node != NULL; j += 1) {
            M[i][j] = cord_node->value;
            cord_node = cord_node->next;
        }
        row_vec = row_vec->next;
    }
    *out_n = n;
    *out_m = m;
    return M;
}

/* Convert dense M[n][m] to array of n cord-rows. */
static void dense_to_cords_rows(double **M, int n, int m, struct cord **out_rows) {
    int i;
    for (i = 0; i < n; i += 1) out_rows[i] = cords_from_array(M[i], m);
}

/* Pack an array of rows into a vector list. */
static struct vector *rows_to_vector(struct cord **rows, int n) {
    int i;
    struct vector *head, *curr, *tail;

    head = NULL; curr = NULL; tail = NULL;

    for (i = 0; i < n; i += 1) {
        tail = (struct vector *)malloc(sizeof(struct vector));
        if (tail == NULL) {
            die();
        }
        tail->cords = rows[i];
        tail->next  = NULL;

        if (head == NULL) {
            head = tail;
            curr = tail;
        } else {
            curr->next = tail;
            curr = tail;
            curr->next = NULL;
        }
    }
    return head;
}

/* Free an n-row dense matrix. (exported in .h) */
void free_matrix(double **M, int n) {
    int i;
    if (M == NULL) return;
    for (i = 0; i < n; i += 1) free(M[i]);
    free(M);
}

/* Free a cord linked list. */
void free_cords(struct cord *head) {
    struct cord *tmp;
    while (head != NULL) {
        tmp = head;
        head = head->next;
        free(tmp);
    }
}

/* Free a vector linked list (including cords). */
void free_vectors(struct vector *head) {
    struct vector *tmp;
    while (head != NULL) {
        tmp = head;
        head = head->next;
        free_cords(tmp->cords);
        free(tmp);
    }
}

/* Return 1 if path ends with ".txt" (case-sensitive; only ".txt" accepted). */
static int has_txt_extension(const char *path) {
    const char *dot = strrchr(path, '.');
    if (!dot) return 0;          /* no dot at all */
    if (dot == path) return 0;   /* filename starts with '.' and nothing before */
    return (dot[1] == 't' && dot[2] == 'x' && dot[3] == 't' && dot[4] == '\0');
}

/* Return 1 if a line has only whitespace. */
static int is_only_ws(const char *buf) {
    const char *p;
    p = buf;
    while (*p != '\0') {
        if (!(*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n')) return 0;
        p += 1;
    }
    return 1;
}

/* Count columns (comma-separated tokens) in a txt line. */
static int count_txt_columns(const char *buf) {
    const char *p;
    int cnt;
    cnt = 1;
    for (p = buf; *p != '\0'; p += 1) if (*p == ',') cnt += 1;
    return cnt;
}

/* Parse a txt line into a new row of size m. */
static int parse_txt_row(const char *line_in, double **out_row, int m) {
    char *line, *tok, *endp;
    int j;
    double v;

    line = my_strdup(line_in);
    if (line == NULL) return -1;

    *out_row = (double *)xmalloc_size((size_t)m, sizeof(**out_row));
    if (*out_row == NULL) { free(line); return -1; }

    tok = strtok(line, ",\r\n");
    j = 0;
    while (tok != NULL && j < m) {
        /* Convert token to double and verify it is fully consumed. */
        v = strtod(tok, &endp);

        /* Skip trailing spaces/tabs after the numeric part. */
        while (*endp == ' ' || *endp == '\t') endp++;

        /* Reject empty tokens or tokens with extra characters. */
        if (endp == tok || *endp != '\0') {
            free(*out_row);
            free(line);
            return -1;
        }

        (*out_row)[j] = v;
        j += 1;
        tok = strtok(NULL, ",\r\n");
    }

    if (j != m) { free(*out_row); free(line); return -1; }
    free(line);
    return 0;
}

/* Grow matrix capacity when needed, keeping M and cap up to date. */
static int grow_matrix(double ***M, int *cap) {
    double **tmp;
    int old;

    old = *cap;

    /* prevent overflow: (*cap * sizeof(*tmp)) * 2 */
    if ((size_t)(*cap) > SIZE_MAX / (2 * sizeof(*tmp))) {
        return -1;
    }

    *cap *= 2;
    tmp = (double **)realloc(*M, (size_t)(*cap) * sizeof(*tmp));
    if (tmp == NULL) { 
        *cap = old; 
        return -1;
    }
    *M = tmp;
    return 0;
}

/* Ensure *buf has at least 'need' bytes (realloc; updates *buf). */
static int ensure_capacity(char **buf, size_t *cap, size_t need) {
    size_t new_cap;
    char *tmp;

    if (need <= *cap) return 0;

    new_cap = (*cap == 0) ? 512 : *cap;
    while (new_cap < need) {
        if (new_cap > SIZE_MAX / 2) { /* doubling would overflow */
            new_cap = need;           /* jump straight to need */
            if (new_cap < *cap) return -1;
            break;
        }
        new_cap *= 2;
    }

    tmp = (char *)realloc(*buf, new_cap);
    if (!tmp) return -1;

    *buf = tmp;
    *cap = new_cap;
    return 0;
}

/* Read a full line from f. */
static char *readline_dynamic(FILE *f) {
    char *buf = NULL;
    size_t cap = 0, len = 0;
    char chunk[4096];

    for (;;) {
        if (fgets(chunk, (int)sizeof(chunk), f) == NULL) {
            if (len == 0) { free(buf); return NULL; }
            if (ensure_capacity(&buf, &cap, len + 1) != 0) { free(buf); return NULL; }
            buf[len] = '\0';
            return buf;
        }

        {
            size_t tlen = strlen(chunk);
            if (ensure_capacity(&buf, &cap, len + tlen + 1) != 0) { free(buf); return NULL; }
            memcpy(buf + len, chunk, tlen);
            len += tlen;
            buf[len] = '\0';
        }

        if (len > 0 && buf[len - 1] == '\n') return buf;
    }
}

/* Read txt into matrix M, setting n and m. 0 on success, -1 on failure. */
static int read_txt(const char *path, double ***out_M, int *out_n, int *out_m) {
    FILE *f;
    int cap, n, m, cols;
    double **M;
    char *line;
    if (!has_txt_extension(path)) return -1;
    cap = 128; n = 0; m = -1;

    f = fopen(path, "r");
    if (f == NULL) return -1;

    M = (double **)xmalloc_size((size_t)cap, sizeof(*M));
    if (M == NULL) { fclose(f); return -1; }

    while ((line = readline_dynamic(f)) != NULL) {
        if (is_only_ws(line)) {
             free(line); free_matrix(M, n); fclose(f); return -1; }
        
        cols = count_txt_columns(line);
        if (m < 0) {
            m = cols;
        } else if (cols != m) {
            free(line); free_matrix(M, n); fclose(f); return -1;}
        
        if (parse_txt_row(line, &M[n], m) != 0) {
            free(line); free_matrix(M, n); fclose(f); return -1;
        }
        free(line);
        n += 1;
        if (n == cap && grow_matrix(&M, &cap) != 0) {
            free_matrix(M, n); fclose(f); return -1;
        }
    }
    fclose(f);
    if (m <= 0 || n <= 0) { free_matrix(M, n); return -1; }
    *out_M = M;
    *out_n = n;
    *out_m = m;
    return 0;
}

/* Print an array of n rows (cord lists), each of length m, with 4 decimals. */
static void print_rows(struct cord **rows, int n, int m) {
    int i, j;
    struct cord *c;

    i = 0;
    j = 0;
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

/* Allocate n×n matrix; clear=1 uses calloc, else malloc. */
static double **alloc_square(int n, int clear) {
    int i ,r;
    double **M;

    M = (double **)xmalloc_size((size_t)n, sizeof(*M));
    if (M == NULL) return NULL;

    for (i = 0; i < n; i += 1) {
        M[i] = clear ? (double *)xcalloc_size((size_t)n, sizeof(**M))
                     : (double *)xmalloc_size((size_t)n, sizeof(**M));
        if (M[i] == NULL) {
            for (r = 0; r < i; r += 1) free(M[r]);
            free(M);
            return NULL;
        }
    }
    return M;
}

/* Compute A[i][*] (j>i) using Gaussian kernel exp(-||xi-xj||^2/2). */
static void build_similarity_row(double **A, double **X, int n, int m, int i) {
    int j, r;
    double s, df, v;
    for (j = i + 1; j < n; j += 1) {
        s = 0.0;
        for (r = 0; r < m; r += 1) {
            df = X[i][r] - X[j][r];
            s += df * df;
        }
        v = exp(-s / 2.0);
        A[i][j] = v;
        A[j][i] = v;
    }
}

/* Build similarity matrix A from dense X. */
static double **build_similarity_matrix(double **X, int n, int m) {
    int i;
    double **A;
    A = alloc_square(n, 1);
    if (A == NULL) return NULL;
    for (i = 0; i < n; i += 1) build_similarity_row(A, X, n, m, i);
    return A;
}

/* Exported: compute A(X) rows. */
void symnmf_sym(struct cord **out_rows, struct vector *data_rows) {
    int n, m;
    double **Xmat, **A;
    Xmat = vectors_to_dense(data_rows, &n, &m);
    if (n == 0 || m == 0) return;
    A = build_similarity_matrix(Xmat, n, m);
    if (A == NULL) { free_matrix(Xmat, n); return; }
    dense_to_cords_rows(A, n, n, out_rows);
    free_matrix(Xmat, n);
    free_matrix(A, n);
}

/* Compute degree array deg[i] = sum_j A[i][j]. */
static double *degrees_from_A(double **A, int n) {
    int i, j;
    double *deg;
    double s;
    deg = (double *)xcalloc_size((size_t)n, sizeof(*deg));
    if (deg == NULL) return NULL;
    for (i = 0; i < n; i += 1) {
        s = 0.0;
        for (j = 0; j < n; j += 1) s += A[i][j];
        deg[i] = s;
    }
    return deg;
}

/* Build diagonal degree matrix D from deg. */
static double **build_D_from_deg(const double *deg, int n) {
    int i;
    double **D;
    D = alloc_square(n, 1);
    if (D == NULL) return NULL;
    for (i = 0; i < n; i += 1) D[i][i] = deg[i];
    return D;
}

/* Helper: build A(X) dense from vectors and return n via out_n. */
static double **build_A_from_vectors(struct vector *data_rows, int *out_n) {
    int n, m, n2;
    struct cord **Arows;
    struct vector *Avec;
    double **Amat;

    n = vectors_len(data_rows);
    Arows = (struct cord **)xcalloc_size((size_t)n, sizeof(*Arows));
    if (Arows == NULL) return NULL;
    symnmf_sym(Arows, data_rows);
    Avec = rows_to_vector(Arows, n);
    Amat = vectors_to_dense(Avec, &n2, &m);
    if (n2 == 0 || n2 != m) { free_vectors(Avec); free(Arows); free_matrix(Amat, n2); return NULL; }
    *out_n = n2;
    free_vectors(Avec);
    free(Arows);
    return Amat;
}

/* Exported: build D internally from X. */
void symnmf_ddg(struct cord **out_rows, struct vector *data_rows) {
    int n;
    double **Amat,**D;
    double *deg;
    
    Amat = build_A_from_vectors(data_rows, &n);
    if (Amat == NULL) return;
    deg = degrees_from_A(Amat, n);
    if (deg == NULL) { free_matrix(Amat, n); return; }
    D = build_D_from_deg(deg, n);
    if (D == NULL) { free(deg); free_matrix(Amat, n); return; }
    dense_to_cords_rows(D, n, n, out_rows);
    free(deg);
    free_matrix(Amat, n);
    free_matrix(D, n);
}

/* Build W = D^{-1/2} A D^{-1/2} from A and deg. */
static double **build_W_from_A_deg(double **A, const double *deg, int n) {
    int i, j;
    double **W;
    double di, dj;
    W = alloc_square(n, 0);   
    if (W == NULL) return NULL;
    for (i = 0; i < n; i += 1) {
        di = (deg[i] > 0.0) ? 1.0 / sqrt(deg[i]) : 0.0;
        for (j = 0; j < n; j += 1) {
            dj = (deg[j] > 0.0) ? 1.0 / sqrt(deg[j]) : 0.0;
            W[i][j] = di * A[i][j] * dj;
        }
    }
    return W;
}

/* Exported: build normalized similarity W internally from X. */
void symnmf_norm(struct cord **out_rows, struct vector *data_rows) {
    int n;
    double **Amat, **W;
    double *deg;
    
    Amat = build_A_from_vectors(data_rows, &n);
    if (Amat == NULL) return;
    deg = degrees_from_A(Amat, n);
    if (deg == NULL) { free_matrix(Amat, n); return; }
    W = build_W_from_A_deg(Amat, deg, n);
    if (W == NULL) { free(deg); free_matrix(Amat, n); return; }
    dense_to_cords_rows(W, n, n, out_rows);
    free(deg);
    free_matrix(Amat, n);
    free_matrix(W, n);
}

/* ========================= SymNMF (multiplicative) ========================= */

/* WH = W * H (sizes n×n and n×k). */
static void W_H(double **W, double **H, double **WH, int n, int k) {
    int i, j, r;
    double Wij;

    for (i = 0; i < n; i += 1) {
        for (r = 0; r < k; r += 1) {
            WH[i][r] = 0.0;
        }

        for (j = 0; j < n; j += 1) {
            Wij = W[i][j];

            if (Wij != 0.0) {
                for (r = 0; r < k; r += 1) {
                    WH[i][r] += Wij * H[j][r];
                }
            }
        }
    }
}

/* HHt = H * H^T. */
static void H_H_T(double **H, double **HHt, int n, int k) {
    int ii, jj, r;
    double s;

    for (ii = 0; ii < n; ii += 1) {
        for (jj = 0; jj < n; jj += 1) {
            s = 0.0;

            for (r = 0; r < k; r += 1) {
                s += H[ii][r] * H[jj][r];
            }

            HHt[ii][jj] = s;
        }
    }
}

/* HHTH = (HHt) * H. */
static void H_H_T_H(double **HHt, double **H, double **HHTH, int n, int k) {
    int i, j, r;
    double s;

    for (i = 0; i < n; i += 1) {
        for (r = 0; r < k; r += 1) {
            s = 0.0;

            for (j = 0; j < n; j += 1) {
                s += HHt[i][j] * H[j][r];
            }

            HHTH[i][r] = s;
        }
    }
}

/* One multiplicative H update step; returns squared Frobenius change. */
static double update_H_step(double **H, double **WH, double **HHTH, int n, int k) {
    int i, r;
    double denom, oldv, newv, diff_sum;
    const double beta = 0.5;
    diff_sum = 0.0;
    for (i = 0; i < n; i += 1) {
        for (r = 0; r < k; r += 1) {
            denom = HHTH[i][r];
            if (denom < TINY_POSITIVE) denom = TINY_POSITIVE;
            oldv = H[i][r];
            newv = oldv * ((1.0 - beta) + beta * (WH[i][r] / denom));
            if (newv < 0.0) newv = 0.0;
            H[i][r] = newv;
            diff_sum += (newv - oldv) * (newv - oldv);
        }
    }
    return diff_sum;
}

/* Copy H cord rows (length k) into dense H[n][k]. */
static double **cords_H_to_dense(struct cord **H, int n, int k) {
    int i;
    int r;
    struct cord *node;
    double **Hmat;
    Hmat = alloc_dense(n, k);
    if (Hmat == NULL) return NULL;
    for (i = 0; i < n; i += 1) {
        node = H[i];
        for (r = 0; r < k && node != NULL; r += 1) {
            Hmat[i][r] = (node->value >= 0.0) ? node->value : 0.0;
            node = node->next;
        }
    }
    return Hmat;
}

/* Allocate SymNMF intermediates (WH, HHt, HHTH). */
static int alloc_intermediates(int n, int k, double ***WH, double ***HHt, double ***HHTH) {
    /* Allocate intermediate matrices:
       - WH    = W * H        (size n × k)   -> numerator in the update
       - HHt   = H * H^T      (size n × n)   -> Gram matrix of H
       - HHTH  = (H H^T) * H  (size n × k)   -> denominator in the update */
    *WH = NULL; *HHt = NULL; *HHTH = NULL;
    
    *WH  = alloc_dense(n, k);
    if (*WH == NULL) return -1;

    *HHt = alloc_square(n, 1);
    if (*HHt == NULL) { free_matrix(*WH, n); *WH = NULL; return -1; }

    *HHTH = alloc_dense(n, k);
    if (*HHTH == NULL) {
        free_matrix(*WH,  n); *WH  = NULL;
        free_matrix(*HHt, n); *HHt = NULL;
        return -1;
    }
    return 0;
}

/* Free SymNMF intermediates. */
static void free_intermediates(double **WH, double **HHt, double **HHTH, int n) {
    free_matrix(WH, n);
    free_matrix(HHt, n);
    free_matrix(HHTH, n);
}

/* Run SymNMF iterations until convergence or MAX_ITER. */
static void run_symnmf_iterations(double **Wmat, double **Hmat,
                                  double **WH, double **HHt, double **HHTH,
                                  int n, int k) {
    int iter;
    double diff_sum;
    /* Main loop */
    for (iter = 0; iter < MAX_ITER; iter += 1) {
        W_H(Wmat, Hmat, WH, n, k);
        H_H_T(Hmat, HHt, n, k);
        H_H_T_H(HHt, Hmat, HHTH, n, k);
        diff_sum = update_H_step(Hmat, WH, HHTH, n, k);
        if (diff_sum < EPS) break;
    }
}

/* Exported: SymNMF multiplicative with beta=0.5; updates H in-place. */
void symnmf_symnmf(struct cord **H, struct vector *W) {
    int n_w, m_w, n ,i, k;
    double **Wmat, **Hmat, **WH, **HHt, **HHTH;

    if (H == NULL || H[0] == NULL || W == NULL) return;
    k = cords_len(H[0]);
    if (k <= 0) return;

    Wmat = vectors_to_dense(W, &n_w, &m_w);
    if (n_w == 0 || n_w != m_w) { free_matrix(Wmat, n_w); return; }
    n = n_w;
    Hmat = cords_H_to_dense(H, n, k);
    if (Hmat == NULL) { free_matrix(Wmat, n); return; }
    if (alloc_intermediates(n, k, &WH, &HHt, &HHTH) != 0) {
        free_matrix(Hmat, n);
        free_matrix(Wmat, n);
        return;
    }
    run_symnmf_iterations(Wmat, Hmat, WH, HHt, HHTH, n, k);
   
    for (i = 0; i < n; i += 1) {
        free_cords(H[i]);
        H[i] = NULL;
    }  
    dense_to_cords_rows(Hmat, n, k, H); 
    
    free_intermediates(WH, HHt, HHTH, n);
    free_matrix(Hmat, n);
    free_matrix(Wmat, n);
}

/* Load input vectors from txt path. */
static int load_input(const char *path, struct vector **out_vecs, double ***out_dense,
                      int *out_n, int *out_d) {
    double **X_dense;
    int n, d;
    struct vector *X_vecs;

    if (read_txt(path, &X_dense, &n, &d) != 0) return -1;

    X_vecs = dense_to_vectors(X_dense, n, d);
    *out_vecs = X_vecs;
    *out_dense = X_dense;
    *out_n = n;
    *out_d = d;
    return 0;
}

/* Run the requested goal and print the resulting matrix. */
static void run_goal_and_print(const char *goal, struct vector *X_vecs, int n) {
    int i;
    struct cord **rows;
    if (strcmp(goal, "symnmf") == 0) die();
    rows = (struct cord **)xcalloc_size((size_t)n, sizeof(*rows));
    if (rows == NULL) die();
    if (strcmp(goal, "sym") == 0) symnmf_sym(rows, X_vecs);
    else if (strcmp(goal, "ddg") == 0) symnmf_ddg(rows, X_vecs);
    else if (strcmp(goal, "norm") == 0) symnmf_norm(rows, X_vecs);
    else { free(rows); die(); }
    print_rows(rows, n, n);
    for (i = 0; i < n; i += 1) free_cords(rows[i]);
    free(rows);
}

/* Simple CLI wrapper. Usage: ./symnmf <goal> <path> */
int main(int argc, char **argv) {
    const char *goal, *path;
    int n, d;
    struct vector *X_vecs;
    double **X_dense;

    if (argc != 3) die();
    goal = argv[1];
    path = argv[2];

    X_vecs = NULL;
    X_dense = NULL;

    if (load_input(path, &X_vecs, &X_dense, &n, &d) != 0) die();
    run_goal_and_print(goal, X_vecs, n);

    free_vectors(X_vecs);
    if (X_dense != NULL) free_matrix(X_dense, n);
    return 0;
}
