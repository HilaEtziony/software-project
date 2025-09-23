#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "symnmf.h"

/*
 * Builds a vector, represented by a C linked list of cords (elements in vector), from a Python list.
 * Returns a pointer to the head of the linked list - first cord in the vector.
 */
static struct cord * build_vector_from_list(PyObject *py_list) {
    int d, i;
    struct cord *head = NULL, *curr = NULL, *new_cord = NULL;
    PyObject* item = NULL;
    double val;
    d = PyObject_Length(py_list);
    if (d <= 0) return NULL;
    for (i = 0; i < d; i++) {
        item = PyList_GetItem(py_list, i);
        if (!item) {  
            free_cords(head);
            return NULL;
        }
        val = PyFloat_AsDouble(item);
        if (PyErr_Occurred()) {  
            free_cords(head);
            PyErr_Clear();
            return NULL;
        }
        new_cord = (struct cord*)malloc(sizeof(struct cord));
        if (!new_cord) {
            free_cords(head);
            return NULL;
        }
        new_cord->value = val;
        new_cord->next = NULL;
        if (head == NULL) {
            head = new_cord;
            curr = head;
        } else {
            curr->next = new_cord;
            curr = curr->next;
        }
    }
    return head;
}

/*
 * Builds a linked list of vectors, from a Python list of vectors (each vector is a Python list).
 * Returns a pointer to the head of the linked list - first vector in the vector's list.
 */
static struct vector * build_vectors_linked_list(PyObject *py_vectors) {
    int n, i;
    struct vector *head = NULL, *curr = NULL, *new_vec = NULL;
    PyObject* py_vec = NULL;

    n = PyObject_Length(py_vectors);
    if (n <= 0) return NULL;
    for (i = 0; i < n; i++) {
        py_vec = PyList_GetItem(py_vectors, i);
        if(!py_vec) {
            free_vectors(head);
            return NULL;
        }
        new_vec = (struct vector*)malloc(sizeof(struct vector));
        if (!new_vec) {
            free_vectors(head);
            return NULL;
        }
        new_vec->cords = build_vector_from_list(py_vec);
        if(!new_vec->cords){
            free_vectors(head);
            return NULL;
        }
        new_vec->next = NULL;
        if (head == NULL) {
            head = new_vec;
            curr = head;
        } else {
            curr->next = new_vec;
            curr = curr->next;
        }
    }
    return head;
}

/*
 * Builds a Python list which represents a vector, from linked list of cords (elements in vector).
 * Returns a pointer to the created Python list.
 */
static PyObject *build_pylist_from_vector(struct cord *cords, int d) {
    int i;
    PyObject *list = NULL, *python_float = NULL;

    list = PyList_New(d);
    if (!list) return NULL;
    for (i = 0; i < d; ++i){
        python_float = PyFloat_FromDouble(cords->value);
        if (!python_float) {
            Py_DECREF(list);
            return NULL;
        }
        if (PyList_SetItem(list, i, python_float) < 0) {
            Py_DECREF(python_float);
            Py_DECREF(list);
            return NULL;
        }
        cords = cords->next;
    }
    return list;
}

/*
 * Convert C matrix to Python list of lists.
 * Returns a pointer to the created Python list of lists.
 */
static PyObject *matrix_to_pylist(struct cord **matrix, int n, int m) {
    PyObject *item = NULL, *result = NULL;
    int i;

    result = PyList_New(n);
    if (!result) return NULL;
    for (i = 0; i < n; i++) {
        item = build_pylist_from_vector(matrix[i], m);
        if (!item) {
            Py_DECREF(result);
            return NULL;
        }
        if (PyList_SetItem(result, i, item) < 0) {
            Py_DECREF(item);
            Py_DECREF(result);
            return NULL;
        }
    }
    return result;
}

/*
 * Free cords matrix.
 */
static void free_cords_matrix(struct cord **matrix, int n) {
    int i;
    for (i = 0; i < n; i++) {
        free_cords(matrix[i]);
    }
    free(matrix);
}

/*
 * Free cords matrix and vectors memory
 */
static void free_cords_matrix_and_vectors(struct cord **matrix, int n, struct vector *vectors) {
    free_cords_matrix(matrix, n);
    free_vectors(vectors);
}

/*
 * Generic goal function to call a C implementation that takes a matrix and a set of vectors as input
 * and produces a matrix accordings to the goal as output.
 * Returns a pointer to the created Python list of lists.
 */
static PyObject *goal(PyObject *args, void (*c_func_goal)(struct cord **, struct vector *)) {
    PyObject *py_X_datapoints = NULL, *result = NULL;
    int n;
    struct vector *X_datapoints = NULL;
    struct cord **matrix = NULL;

    /* Parse Args */
    if (!PyArg_ParseTuple(args, "O", &py_X_datapoints)) return NULL;
    n = PyObject_Length(py_X_datapoints);
    if (n <= 0) return NULL;

    /* Process X datapoint */
    X_datapoints = build_vectors_linked_list(py_X_datapoints);
    if (!X_datapoints) return NULL;

    matrix = (struct cord**)calloc(n, sizeof(struct cord*));
    if (!matrix) {
        free_vectors(X_datapoints);
        return NULL;
    }

    /* Calling the C implementation */
    c_func_goal(matrix, X_datapoints);
    /* Conversion matrix from C to Python (there is no need to check result because the memory cleanup happens next anyway) */
    result = matrix_to_pylist(matrix, n, n);
    /* Memory cleanup */
    free_cords_matrix_and_vectors(matrix, n, X_datapoints);

    return result;
}

/*
 * Build an array of vectors from a Python list of lists.
 * Returns a pointer to the created array of vectors.
 */
static struct cord ** build_cords_matrix_from_pylist(PyObject *py_matrix, int n) {
    PyObject *py_cords_list = NULL;
    int i;

    struct cord **cords_matrix = (struct cord **)calloc(n, sizeof(struct cord *));
    if (!cords_matrix) return NULL;

    for (i = 0; i < n; i++) {
        py_cords_list = PyList_GetItem(py_matrix, i);
        if(py_cords_list == NULL) {
            free_cords_matrix(cords_matrix, i);
            return NULL;
        }
        cords_matrix[i] = build_vector_from_list(py_cords_list);
        if (cords_matrix[i] == NULL) {
            free_cords_matrix(cords_matrix, i);
            return NULL;
        }
    }

    return cords_matrix;
}

/* 
 * Run symNMF full algorithm from C and return the resulting H matrix.
 */
static PyObject *symnmf(PyObject *Py_UNUSED(self), PyObject *args) {
    PyObject *py_H_matrix = NULL, *py_W_matrix = NULL, *result = NULL;
    int k, n;
    struct vector *W_matrix = NULL;
    struct cord **H_matrix = NULL;

    if (!PyArg_ParseTuple(args, "OO", &py_H_matrix, &py_W_matrix)) return NULL;
    n = PyObject_Length(py_H_matrix);
    if (n <= 0) return NULL;

    W_matrix = build_vectors_linked_list(py_W_matrix);
    if (!W_matrix) return NULL;
    H_matrix = build_cords_matrix_from_pylist(py_H_matrix, n);
    if (!H_matrix) {
        free_vectors(W_matrix);
        return NULL;
    }

    /* Calling the C implementation */
    symnmf_symnmf(H_matrix, W_matrix);
    /* Get number of cords of a vector in H*/
    k = cords_len(H_matrix[0]);
    /* Conversion matrix from C to Python (there is no need to check result because the memory cleanup happens next anyway) */
    result = matrix_to_pylist(H_matrix, n, k);
    /* Memory cleanup */
    free_cords_matrix_and_vectors(H_matrix, n, W_matrix);

    return result;
}

/*
 * Calculate and return the similarity matrix.
 */
static PyObject * sym(PyObject *Py_UNUSED(self), PyObject *args) {
    return goal(args, symnmf_sym);
}

/*
 * Calculate and return the diagonal degree matrix.
 */
static PyObject * ddg(PyObject *Py_UNUSED(self), PyObject *args) {
    return goal(args, symnmf_ddg);
}

/*
 * Calculate and return the normalized similarity matrix.
 */
static PyObject * norm(PyObject *Py_UNUSED(self), PyObject *args) {
    return goal(args, symnmf_norm);
}

/* === Module Method Table === */
static PyMethodDef symnmfMethods[] = {
    {
        "symnmf",  
        (PyCFunction)symnmf,          
        METH_VARARGS,
        PyDoc_STR("Run symNMF full algorithm from C")
    },{
        "sym",  
        (PyCFunction)sym,           
        METH_VARARGS,
        PyDoc_STR("Calculate and return the similarity matrix")
    },{
        "ddg",  
        (PyCFunction)ddg,  
        METH_VARARGS,
        PyDoc_STR("Calculate and return the diagonal degree matrix")
    },{
        "norm",  
        (PyCFunction)norm,           
        METH_VARARGS,
        PyDoc_STR("Calculate and return the normalized similarity matrix")
    },{
        NULL, NULL, 0, NULL
    }
};

/* === Module Definition === */
static struct PyModuleDef symnmfmodule = {
    PyModuleDef_HEAD_INIT,
    "symnmfmodule", /* name of module */
    NULL,
    -1,
    symnmfMethods, 
    NULL,
    NULL,
    NULL,
    NULL
};

/* === Module Initialization === */
PyMODINIT_FUNC PyInit_symnmfmodule(void) {
    PyObject *projectmodule;
    projectmodule = PyModule_Create(&symnmfmodule);
    if (!projectmodule) {
        return NULL;
    }
    return projectmodule;
}
