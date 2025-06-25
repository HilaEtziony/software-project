# define PY_SSIZE_T_CLEAN
# include <Python.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "kmeans.h"



static struct cord* build_cord_from_list(PyObject* py_list) {
    int d, i;
    struct cord* head = NULL;
    struct cord* curr = NULL;
    struct cord* new_cord = NULL;
    PyObject* item = NULL;
    double val;

    d = PyObject_Length(py_list);
    for (i = 0; i < d; i++) {
        item = PyList_GetItem(py_list, i);
        val = PyFloat_AsDouble(item);
        new_cord = (struct cord*)malloc(sizeof(struct cord));
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

/* Constructs a C linked list as taught in class from Python multidimensional array */
static struct vector* build_vector_list(PyObject* py_vectors) {
    int n, i;
    struct vector* head = NULL;
    struct vector* curr = NULL;
    struct vector* new_vec = NULL;
    PyObject* py_vec = NULL;

    n = PyObject_Length(py_vectors);
    for (i = 0; i < n; i++) {
        py_vec = PyList_GetItem(py_vectors, i);
        new_vec = (struct vector*)malloc(sizeof(struct vector));
        new_vec->cords = build_cord_from_list(py_vec);
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

static PyObject* cords_to_pylist(struct cord* c , int d) {
    int i;
    PyObject* list = NULL;
    PyObject* python_float = NULL;

    list = PyList_New(d);
    for (i = 0; i < d; ++i){
        python_float = PyFloat_FromDouble(c->value);
        PyList_SetItem(list, i, python_float);
        c = c->next;
    }
    return list;
}

/* Wrapping the function */

static PyObject* fit(PyObject* self, PyObject* args) {
    PyObject *py_centroids = NULL, *py_vectors = NULL;
    PyObject* py_c = NULL;
    PyObject* c = NULL;
    PyObject* result = NULL;
    PyObject* first = NULL;

    int k, iter, i, d;
    double epsilon;

    struct vector* vec_list = NULL;
    struct cord** centroids = NULL;

    if (!PyArg_ParseTuple(args, "idOO", &iter, &epsilon, &py_centroids, &py_vectors)) {
        return NULL;
    }

    k = PyObject_Length(py_centroids);
    if (k < 0) {
        return NULL;
    }

    first = PyList_GetItem(py_centroids, 0);
    d = PyObject_Length(first);
    if (d < 0) {
        return NULL;
    }

    vec_list = build_vector_list(py_vectors);

    centroids = (struct cord**)calloc(k, sizeof(struct cord*));
    for (i = 0; i < k; i++) {
        py_c = PyList_GetItem(py_centroids, i);
        centroids[i] = build_cord_from_list(py_c);
    }

    /* Calling the C implementation */
    kmeans_fit(k, iter, epsilon, centroids, vec_list);
    
    /* New conversion from C to Python */
    result = PyList_New(k);
    for (i = 0; i < k; i++) {
        c = cords_to_pylist(centroids[i], d);
        PyList_SetItem(result, i, c); 
    }

    /* Memory cleanup */
    for (i = 0; i < k; i++) {
        free_cords(centroids[i]);
    }
    free(centroids);
    free_vectors(vec_list);

    return result;
}

/* 3. === Module Method Table === */

static PyMethodDef kmeansMethods[] = {
    {"fit",  /* the Python method name that will be used */
        (PyCFunction)fit,  /* the C-function that implements the Python function and returns static PyObject*  */         
        METH_VARARGS,/* flags indicating parameters
accepted for this function */
        PyDoc_STR("Run KMeans algorithm from C")},
    {NULL, NULL, 0, NULL}
};

/* 4. === Module Definition === */

static struct PyModuleDef kmeansmodule = {
    PyModuleDef_HEAD_INIT,
    "mykmeanspp", /* name of module */
    NULL,
    -1,
    kmeansMethods /* the PyMethodDef array from before 
    containing the methods of the extension */
};

/* 5. === Module Initialization === */

PyMODINIT_FUNC PyInit_mykmeanspp(void) {
    PyObject *m;
    m = PyModule_Create(&kmeansmodule);
    if (!m) {
        return NULL;
    }
    return m;
}
