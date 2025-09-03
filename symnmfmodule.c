# define PY_SSIZE_T_CLEAN
# include <Python.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "symnmf.h"



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

static PyObject* symnmf(PyObject* self, PyObject* args) {
    PyObject *py_H_metrix = NULL, *py_W_metrix = NULL;
    PyObject* py_c = NULL;
    PyObject* c = NULL;
    PyObject* result = NULL;
    PyObject* first = NULL;
    int k, i, n;
    struct vector* W_metrix_lines_list = NULL;
    struct cord** H_metrix_lines_list = NULL;

    if (!PyArg_ParseTuple(args, "iOO", &k, &py_H_metrix, &py_W_metrix)) {
        return NULL;
    }

    n = PyObject_Length(py_H_metrix);
    if (n < 0) {
        return NULL;
    }

    W_metrix_lines_list = build_vector_list(py_W_metrix);

    H_metrix_lines_list = (struct cord**)calloc(n, sizeof(struct cord*));
    for (i = 0; i < n; i++) {
        py_c = PyList_GetItem(py_H_metrix, i);
        H_metrix_lines_list[i] = build_cord_from_list(py_c);
    }

    /* Calling the C implementation */
    symnmf_algo(k, H_metrix_lines_list, W_metrix_lines_list);
    
    /* New conversion from C to Python */
    result = PyList_New(n);
    for (i = 0; i < n; i++) {
        c = cords_to_pylist(H_metrix_lines_list[i], k);
        PyList_SetItem(result, i, c); 
    }

    /* Memory cleanup */
    for (i = 0; i < k; i++) {
        free_cords(H_metrix_lines_list[i]);
    }
    free(H_metrix_lines_list);
    free_vectors(H_metrix_lines_list);

    return result;
}

static PyObject* sym(PyObject* self, PyObject* args) {
    PyObject *py_X_datapoints = NULL;
    PyObject* py_c = NULL;
    PyObject* c = NULL;
    PyObject* result = NULL;
    PyObject* first = NULL;
    int d, i, n;
    struct vector* W_metrix_lines_list = NULL;
    struct cord** H_metrix_lines_list = NULL;

    if (!PyArg_ParseTuple(args, "O", &py_X_datapoints)) {
        return NULL;
    }

    d = PyObject_Length(py_X_datapoints);
    if (d < 0) {
        return NULL;
    }

    W_metrix_lines_list = build_vector_list(py_X_datapoints);

    H_metrix_lines_list = (struct cord**)calloc(n, sizeof(struct cord*));

    /* Calling the C implementation */
    sym_algo(H_metrix_lines_list, W_metrix_lines_list);
    
    /* New conversion from C to Python */
    result = PyList_New(n);
    for (i = 0; i < n; i++) {
        c = cords_to_pylist(H_metrix_lines_list[i], n);
        PyList_SetItem(result, i, c); 
    }

    /* Memory cleanup */
    for (i = 0; i < n; i++) {
        free_cords(H_metrix_lines_list[i]);
    }
    free(H_metrix_lines_list);
    free_vectors(H_metrix_lines_list);

    return result;
}

//TODO - after gad will write his part, check symnmf and sym are fir and duplicate sym to ddg and norm

/* 3. === Module Method Table === */

static PyMethodDef symnmfMethods[] = {
    {
        "symnmf",  /* the Python method name that will be used */
        (PyCFunction)symnmf,  /* the C-function that implements the Python function and returns static PyObject*  */         
        METH_VARARGS,/* flags indicating parameters accepted for this function */
        PyDoc_STR("Run symNMF full algorithm from C")
    },{
        "sym",  /* the Python method name that will be used */
        (PyCFunction)sym,  /* the C-function that implements the Python function and returns static PyObject*  */         
        METH_VARARGS,/* flags indicating parameters accepted for this function */
        PyDoc_STR("Calculate and output the similarity matrix")
    },{
        "ddg",  /* the Python method name that will be used */
        (PyCFunction)ddg,  /* the C-function that implements the Python function and returns static PyObject*  */         
        METH_VARARGS,/* flags indicating parameters accepted for this function */
        PyDoc_STR("Calculate and output the diagonal degree matrix")
    },{
        "norm",  /* the Python method name that will be used */
        (PyCFunction)norm,  /* the C-function that implements the Python function and returns static PyObject*  */         
        METH_VARARGS,/* flags indicating parameters accepted for this function */
        PyDoc_STR("Calculate and output the normalized similarity matrix")
    },{
        NULL, NULL, 0, NULL
    }
};

/* 4. === Module Definition === */

static struct PyModuleDef symnmfmodule = {
    PyModuleDef_HEAD_INIT,
    "symnmfmodule", /* name of module */
    NULL,
    -1,
    symnmfMethods /* the PyMethodDef array from before containing the methods of the extension */
};

/* 5. === Module Initialization === */

PyMODINIT_FUNC PyInit_symnmfmodule(void) {
    PyObject *projectmodule;
    projectmodule = PyModule_Create(&symnmfmodule);
    if (!projectmodule) {
        return NULL;
    }
    return projectmodule;
}
