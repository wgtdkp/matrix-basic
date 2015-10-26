#include <Python.h>
#include <structmember.h>
#include "matrix.h"
#include "direct.h"
#include "iterative.h"
#include "eigenvalue.h"
#include "norm.h"

#define LIST(L, i)      PyList_GetItem((L), (i))

static const char mym_doc[] = 
    "mym is a simple matrix library, it supports\n"
    "many basical matrix-based operations.";


typedef struct{
    PyObject_HEAD
    Matrix* M;
} MYM_Matrix;

static inline int mym_matrix_converter(PyObject* obj_m, Matrix** mat);
static inline int mym_bool_converter(PyObject* obj_b, bool* b);
static MYM_Matrix* mym_matrix_create(Matrix* M);
static PyObject* mym_matrix_new(PyTypeObject* type, PyObject* args, PyObject *kwds);
static void mym_matrix_dealloc(MYM_Matrix* self);
static int mym_matrix_init(MYM_Matrix* self, PyObject *args, PyObject* kwds);
static PyObject* mym_matrix_size(MYM_Matrix* self);
static PyObject* mym_matrix_copy(MYM_Matrix* self);
static PyObject* mym_matrix_repr(MYM_Matrix* self);
static PyObject* mym_eye(PyObject* self, PyObject* args, PyObject* kwds);
static PyObject* mym_add(PyObject* self, PyObject* args, PyObject* kwds);
static PyObject* mym_sub(PyObject* self, PyObject* args, PyObject* kwds);
static PyObject* mym_mul(PyObject* self, PyObject* args, PyObject* kwds);
static PyObject* mym_mul_cons(PyObject* self, PyObject* args, PyObject* kwds);
static PyObject* mym_det(PyObject* self, PyObject* args, PyObject* kwds);
static PyObject* mym_trans(PyObject* self, PyObject* args, PyObject* kwds);
static PyObject* mym_inv(PyObject* self, PyObject* args, PyObject* kwds);
static PyObject* mym_gauss_elim(PyObject* self, PyObject* args, PyObject* kwds);
static PyObject* mym_me_gauss_elim(PyObject* self, PyObject* args, PyObject* kwds);
static PyObject* mym_tri_decomp(PyObject* self, PyObject* args, PyObject* kwds);
static PyObject* mym_me_tri_decomp(PyObject* self, PyObject* args, PyObject* kwds);
static PyObject* mym_cholesky_decomp(PyObject* self, PyObject* args, PyObject* kwds);
static PyObject* mym_en_cholesky_decomp(PyObject* self, PyObject* args, PyObject* kwds);
static PyObject* mym_chasing_method(PyObject* self, PyObject* args, PyObject* kwds);
static PyObject* mym_jacobi_iter(PyObject* self, PyObject* args, PyObject* kwds);
static PyObject* mym_gauss_seidel_iter(PyObject* self, PyObject* args, PyObject* kwds);
static PyObject* mym_sor(PyObject* self, PyObject* args, PyObject* kwds);
static PyObject* mym_pow_method(PyObject* self, PyObject* args, PyObject* kwds);
static PyObject* mym_norm(PyObject* self, PyObject* args, PyObject* kwds);

PyMODINIT_FUNC PyInit_mym(void);



static PyMemberDef mym_matrix_members[] = {

};

static PyMethodDef mym_matrix_methods[] = {
    {"size", (PyCFunction)mym_matrix_size, METH_VARARGS, PyDoc_STR("get the matrix size")},
    {"copy", (PyCFunction)mym_matrix_copy, METH_VARARGS, PyDoc_STR("make a deep copy")},
    {NULL}
};

static PyTypeObject MYM_Matrix_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "mym.Matrix",             /* tp_name */
    sizeof(MYM_Matrix), /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor)mym_matrix_dealloc, /* tp_dealloc */
    0,                         /* tp_print */
    0,                         /* tp_getattr */
    0,                         /* tp_setattr */
    0,                         /* tp_reserved */
    (reprfunc)mym_matrix_repr,                         /* tp_repr */
    0,                         /* tp_as_number */
    0,                         /* tp_as_sequence */
    0,                         /* tp_as_mapping */
    0,                         /* tp_hash  */
    0,                         /* tp_call */
    0,                         /* tp_str */
    0,                         /* tp_getattro */
    0,                         /* tp_setattro */
    0,                         /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,        /* tp_flags */
    "Matrix objects",           /* tp_doc */
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    mym_matrix_methods,             /* tp_methods */
    mym_matrix_members,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)mym_matrix_init,      /* tp_init */
    0,                         /* tp_alloc */
    mym_matrix_new,                 /* tp_new */
};

static PyMethodDef mym_methods[] = {
    {"eye",     (PyCFunction)mym_eye, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("create eye.")},
    {"add",     (PyCFunction)mym_add, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("add two Matrix.")},
    {"sub",     (PyCFunction)mym_sub, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("make substraction.")},
    {"mul",     (PyCFunction)mym_mul, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("mul two Matrix.")},
    {"mul_cons", (PyCFunction)mym_mul_cons, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("mul a constant.")},
    {"det",     (PyCFunction)mym_det, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("clculate the determinant.")},
    {"trans",   (PyCFunction)mym_trans, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("get the transpose matrix.")},
    {"inv",     (PyCFunction)mym_inv, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("get the inverse matrix.")},
    {"gauss_elim", (PyCFunction)mym_gauss_elim, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("do Guass elimination")},
    {"me_gauss_elim", (PyCFunction)mym_me_gauss_elim, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("do main element gauss elimination")},
    {"tri_decomp",  (PyCFunction)mym_tri_decomp, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("do triangle decomposition")},
    {"me_tri_decomp", (PyCFunction)mym_me_tri_decomp, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("do main elment triangle decomposition")},
    {"cholesky_decomp", (PyCFunction)mym_cholesky_decomp, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("do Cholesky decomposition")},
    {"en_cholesky_decomp", (PyCFunction)mym_en_cholesky_decomp, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("do enhanced Cholesky decomposition")},
    {"chasing_method", (PyCFunction)mym_chasing_method, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("do Chasing method")},
    {"jacobi", (PyCFunction)mym_jacobi_iter, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("do Jacobi iteration")},
    {"gauss_seidel", (PyCFunction)mym_gauss_seidel_iter, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("do Gauss-Seidel iteration")},
    {"sor", (PyCFunction)mym_sor, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("do SOR iteration")},
    {"pow_method", (PyCFunction)mym_pow_method, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("calc biggest eigenvalue with pow method")},
    {"norm", (PyCFunction)mym_norm, METH_VARARGS | METH_KEYWORDS, PyDoc_STR("calc matrix norm")},
    {NULL}
};

static struct PyModuleDef mym_module = {
    PyModuleDef_HEAD_INIT,
    .m_name = "mym",
    .m_doc = mym_doc,
    -1,
    .m_methods = mym_methods,
    NULL, NULL, NULL, NULL
};

#define IS_MATRIX(m)    (PyObject_IsInstance((m), (PyObject*)(&MYM_Matrix_Type)))
#define IS_DOUBLE(x)    (PyFloat_Check((x)) || PyLong_Check((x)))

static inline int mym_matrix_converter(PyObject* obj_m, Matrix** M)
{
    Matrix* tmp;
    if(1 != IS_MATRIX(obj_m)) {
        PyErr_SetString(PyExc_TypeError, "the arg should be mym.Matrix");
        return 0;
    }
    if(IS_MATRIX(obj_m)) {
        tmp = ((MYM_Matrix*)obj_m)->M;
        if(!IS_VECTOR(tmp) && !IS_SQUARE(tmp)) {
            PyErr_SetString(PyExc_ValueError, \
                "only square matrix and vector is supported now");
            return 0;
        }
        *M = tmp;
    }
    return 1;
}

static inline int mym_bool_converter(PyObject* obj_b, bool* b)
{
    if(!PyBool_Check(obj_b))
        return 0;
    *b = (Py_True == obj_b);
    return 1;
}

static MYM_Matrix* mym_matrix_create(Matrix* M)
{
    MYM_Matrix* mat;
    if(NULL == M) {
        PyErr_SetString(PyExc_ValueError, "can't build matrix");
        return NULL;
    }
    mat = (MYM_Matrix*)mym_matrix_new(&MYM_Matrix_Type, NULL, NULL);
    mat->M = M;
    return mat;
}

static PyObject* mym_matrix_new(PyTypeObject* type, PyObject* args, PyObject *kwds)
{
    MYM_Matrix* self;
    self = (MYM_Matrix*)type->tp_alloc(type, 0);
    if(NULL != self)
        self->M = NULL;
    return (PyObject*)self;
}



static void mym_matrix_dealloc(MYM_Matrix* self)
{
    free(self->M);
    self->M = NULL;
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static int mym_matrix_init(MYM_Matrix* self, PyObject *args, PyObject* kwds)
{
    static char* kwlist[] = {"val", NULL};
    int m, n, i, j;
    PyObject *val = NULL, *val_item, *data;
    
    //arguments checking
    if(!PyArg_ParseTupleAndKeywords(args, kwds, "O", kwlist, &val))
        goto VALUE_ERROR;
    if(!PyList_Check(val) || PyList_Size(val) < 1)
        goto VALUE_ERROR;
    m = PyList_Size(val);
    val_item = LIST(val, 0);
    if(!PyList_Check(val_item) || PyList_Size(val_item) < 1)
        goto VALUE_ERROR;
    n = PyList_Size(val_item);
    for(i = 0; i < m; i++) {
        val_item = LIST(val, i);
        if(!PyList_Check(val_item) || PyList_Size(val_item) != n)
            goto VALUE_ERROR;
        for(j = 0; j < n; j++) {
            data = LIST(val_item, j);
            if(!IS_DOUBLE(data))
                goto VALUE_ERROR;
        }
    }

    //init the Matrix
    self->M = create_matrix(m, n);
    for(i = 0; i < m; i++)
        for(j = 0; j < n; j++){
            data = LIST(LIST(val, i), j);
            if(PyFloat_Check(data))
                self->M->mem[i][j] = PyFloat_AsDouble(data);
            else
                self->M->mem[i][j] = (double)PyLong_AsLong(data); 
        }
    return 0;

VALUE_ERROR:
    PyErr_SetString(PyExc_ValueError, "the 'val' is bad type or format");
     return -2;
}

static PyObject* mym_matrix_size(MYM_Matrix* self)
{
    Matrix* M = self->M;
    if(NULL == M) {
        PyErr_SetString(PyExc_AttributeError, "Matrix is empty");
        return NULL;
    }
    return Py_BuildValue("[ii]", M->m, M->n);
}

static PyObject* mym_matrix_copy(MYM_Matrix* self)
{
    Matrix *M = self->M, *M_copy = NULL;
    MYM_Matrix* ret = NULL;
    if(NULL == M) {
        PyErr_SetString(PyExc_AttributeError, "Matrix is empty");
        return NULL;
    }
    M_copy = copy(M);
    ret = mym_matrix_create(M_copy);
    return (PyObject*)ret;
}

static PyObject* mym_matrix_repr(MYM_Matrix* self)
{
    Matrix* M = self->M;
    int i, j;
    PyObject* repr = PyUnicode_FromString("");
    static char* buffer = NULL;
    if(NULL == buffer)
        buffer = (char*)malloc(sizeof(char) * 200);
    for(i = 0; i < M->m; i++) {
        PyObject* row = PyUnicode_FromString("[");
        for(j = 0; j < M->n; j++) {
            sprintf(buffer, "%.8lf\t", M->mem[i][j]);
            buffer[199] = 0;
            row = PyUnicode_Concat(row, PyUnicode_FromString(buffer));
        }
        row = PyUnicode_Concat(row, PyUnicode_FromString("]\n"));
        repr = PyUnicode_Concat(repr, row);
    }
    return repr;
}

static PyObject* mym_eye(PyObject* self, PyObject* args, PyObject* kwds)
{
    int n;
    static char* kwlist[] = {"n", NULL};
    if(!PyArg_ParseTupleAndKeywords(args, kwds, "i:eye", kwlist, &n))
        return NULL;
    return (PyObject*)mym_matrix_create(create_eye(n));
 }

#define COMMON_ARGS_PARSING_A(args, kwds, name) Matrix* A;\
    bool step = false;\
    static char* kwlist[] = {"A", "step", NULL};\
    if(!PyArg_ParseTupleAndKeywords(args, kwds, \
        "O&|O&:"name, kwlist, \
        &mym_matrix_converter, &A, \
        &mym_bool_converter, &step)) {\
        return NULL;\
    }

#define COMMON_ARGS_PARSING_A_B(args, kwds, name) Matrix *A, *B;\
    bool step = false;\
    Matrix* M;\
    static char* kwlist[] = {"A", "B", "step", NULL};\
    if(!PyArg_ParseTupleAndKeywords((args), (kwds), \
        "O&O&|O&:"name, kwlist, \
        &mym_matrix_converter, &A, \
        &mym_matrix_converter, &B, \
        &mym_bool_converter, &step)) {\
        return NULL;\
    }

static PyObject* mym_mul(PyObject* self, PyObject* args, PyObject* kwds)
{
    COMMON_ARGS_PARSING_A_B(args, kwds, "mul")
    if(A->n != B->m) {
        PyErr_SetString(PyExc_ValueError, \
            "can't mul the two Matrixes");
        return NULL;
    }

    M = mul(A, B);
    return (PyObject*)mym_matrix_create(M);
}


static PyObject* mym_mul_cons(PyObject* self, PyObject* args, PyObject* kwds)
{
    Matrix* A;
    double x;
    bool step = false;
    static char* kwlist[] = {"A", "x", "step", NULL};
    if(!PyArg_ParseTupleAndKeywords((args), (kwds), \
        "O&d|O&:mul_cons", kwlist, \
        &mym_matrix_converter, &A, \
        &x, \
        &mym_bool_converter, &step)) {
        return NULL;
    }
    return (PyObject*)mym_matrix_create(mul_cons(A, x));
}


static PyObject* mym_add(PyObject* self, PyObject* args, PyObject* kwds)
{
    COMMON_ARGS_PARSING_A_B(args, kwds, "add")
    if(!IS_SAME_SIZE(A, B)) {
        PyErr_SetString(PyExc_ValueError, \
            "can't add the two Matrixes, they have different size");
        return NULL;
    }

    M = add(A, B);
    return (PyObject*)mym_matrix_create(M);
}

static PyObject* mym_sub(PyObject* self, PyObject* args, PyObject* kwds)
{
    COMMON_ARGS_PARSING_A_B(args, kwds, "sub")
    if(!IS_SAME_SIZE(A, B)) {
        PyErr_SetString(PyExc_ValueError, \
            "can't make substraction, they have different size");
        return NULL;
    }

    M = sub(A, B);
    return (PyObject*)mym_matrix_create(M);
}




static PyObject* mym_det(PyObject* self, PyObject* args, PyObject* kwds)
{
    COMMON_ARGS_PARSING_A(args, kwds, "det")
    if(!IS_SQUARE(A)) {
        PyErr_SetString(PyExc_ValueError, \
            "only square matrix has determinant"); 
        return NULL;
    }
    return PyFloat_FromDouble(det(A));
}

static PyObject* mym_trans(PyObject* self, PyObject* args, PyObject* kwds)
{
    COMMON_ARGS_PARSING_A(args, kwds, "trans")
    return (PyObject*)mym_matrix_create(transpose(A));
}

static PyObject* mym_inv(PyObject* self, PyObject* args, PyObject* kwds)
{
    Matrix* inv;
    COMMON_ARGS_PARSING_A(args, kwds, "inv")
    inv = inverse(A);
    if(NULL == inv) {
        PyErr_SetString(PyExc_ValueError, \
            "the matrix is a singular matrix, no inervse matrix");
        return NULL;
    }
    return (PyObject*)mym_matrix_create(inv);
}

    //EQUATION_CHECK(A, B)

#define EQUATION_CHECK(A, B)    if(!IS_EQUATION((A), (B))) {\
        PyErr_SetString(PyExc_ValueError, \
            "the two matrixes is not a legal equation");\
        return NULL;\
    }

static PyObject* mym_gauss_elim(PyObject* self, PyObject* args, PyObject* kwds)
{
    COMMON_ARGS_PARSING_A_B(args, kwds, "gauss_elim")
    EQUATION_CHECK(A, B)
    M = gauss_elim(A, B, step);
    return (PyObject*)mym_matrix_create(M);
}

static PyObject* mym_me_gauss_elim(PyObject* self, PyObject* args, PyObject* kwds)
{
    COMMON_ARGS_PARSING_A_B(args, kwds, "me_gauss_elim")
    EQUATION_CHECK(A, B)
    M = me_gauss_elim(A, B, step);
    return (PyObject*)mym_matrix_create(M);
}

static PyObject* mym_tri_decomp(PyObject* self, PyObject* args, PyObject* kwds)
{
    COMMON_ARGS_PARSING_A_B(args, kwds, "tri_decomp")
    EQUATION_CHECK(A, B)
    if(!is_ordered_main_subdet(A, &ne, 0)) {
        PyErr_SetString(PyExc_ValueError, \
            "one of A matrix's ordered main sub determinant is zero, "
            "thus A can't be decomposed");
        return NULL;
    }
    M = tri_decomp(A, B, step);
    return (PyObject*)mym_matrix_create(M);
}

static PyObject* mym_me_tri_decomp(PyObject* self, PyObject* args, PyObject* kwds)
{
    COMMON_ARGS_PARSING_A_B(args, kwds, "me_tri_decomp")
    EQUATION_CHECK(A, B)
    if(!is_ordered_main_subdet(A, &ne, 0)) {
        PyErr_SetString(PyExc_ValueError, \
            "one of A matrix's ordered main sub determinant is zero, "
            "thus A can't be decomposed");
        return NULL;
    }
    M = me_tri_decomp(A, B, step);
    return (PyObject*)mym_matrix_create(M);
}

static PyObject* mym_cholesky_decomp(PyObject* self, PyObject* args, PyObject* kwds)
{
    COMMON_ARGS_PARSING_A_B(args, kwds, "cholesky_decomp")
    EQUATION_CHECK(A, B)
    if(!(is_symmetrical(A) && is_ordered_main_subdet(A, &gt, 0))) {
        PyErr_SetString(PyExc_ValueError, \
            "the A matrix is not a symmetric positive definite");
        return NULL;
    }
    M = cholesky_decomp(A, B, step);
    return (PyObject*)mym_matrix_create(M);
}

static PyObject* mym_en_cholesky_decomp(PyObject* self, PyObject* args, PyObject* kwds)
{
    COMMON_ARGS_PARSING_A_B(args, kwds, "en_cholesky_decomp")
    EQUATION_CHECK(A, B)
    if(!(is_symmetrical(A) && is_ordered_main_subdet(A, &gt, 0))) {
        PyErr_SetString(PyExc_ValueError, \
            "the A matrix is not a symmetric positive definite");
        return NULL;
    }
    M = en_cholesky_decomp(A, B, step);
    return (PyObject*)mym_matrix_create(M);
}

static PyObject* mym_chasing_method(PyObject* self, PyObject* args, PyObject* kwds)
{
    COMMON_ARGS_PARSING_A_B(args, kwds, "chasing_method")
    EQUATION_CHECK(A, B)
    if(!is_diagnoal_dominance_3(A)) {
        PyErr_SetString(PyExc_ValueError, \
            "the matrix is not a 3 diagnoal dominance");
        return NULL;
    }
    M = chasing_method(A, B, step);
    return (PyObject*)mym_matrix_create(M);
}

static PyObject* mym_jacobi_iter(PyObject* self, PyObject* args, PyObject* kwds)
{
    Matrix *A, *B;
    bool step = false;
    Matrix* M;
    double delta;
    static char* kwlist[] = {"A", "B", "delta", "step", NULL};
    if(!PyArg_ParseTupleAndKeywords(args, kwds, 
        "O&O&d|O&:jacobi", kwlist, 
        &mym_matrix_converter, &A, 
        &mym_matrix_converter, &B,
        &delta, &mym_bool_converter, &step)) {
        return NULL;
    }
    EQUATION_CHECK(A, B)
    if(!is_jacobi_convergent(A, step)) {
        PyErr_SetString(PyExc_ValueError, 
            "the jacobi iteration is not convergent, "
            "set 'step' True to get more information");
        return NULL;
    }
    M = jacobi_iter(A, B, delta, step);
    return (PyObject*)mym_matrix_create(M);
}

static PyObject* mym_gauss_seidel_iter(PyObject* self, PyObject* args, PyObject* kwds)
{
    Matrix *A, *B;
    bool step = false;
    Matrix* M;
    double delta;
    static char* kwlist[] = {"A", "B", "delta", "step", NULL};
    if(!PyArg_ParseTupleAndKeywords(args, kwds, 
        "O&O&d|O&:gauss_seidel", kwlist, 
        &mym_matrix_converter, &A, 
        &mym_matrix_converter, &B,
        &delta, &mym_bool_converter, &step)) {
        return NULL;
    }
    EQUATION_CHECK(A, B)
    if(!is_gauss_seidel_convergent(A, step)) {
        PyErr_SetString(PyExc_ValueError, 
            "the jacobi iteration is not convergent, "
            "set 'step' True to get more information");
        return NULL;
    }
    M = gauss_seidel_iter(A, B, delta, step);
    return (PyObject*)mym_matrix_create(M);
}

static PyObject* mym_sor(PyObject* self, PyObject* args, PyObject* kwds)
{
    Matrix *A, *B;
    bool step = false;
    Matrix* M;
    double w, delta;
    static char* kwlist[] = {"A", "B", "w", "delta", "step", NULL};
    if(!PyArg_ParseTupleAndKeywords(args, kwds, 
        "O&O&dd|O&:sor", kwlist, 
        &mym_matrix_converter, &A, 
        &mym_matrix_converter, &B,
        &w, &delta, &mym_bool_converter, &step)) {
        return NULL;
    }
    EQUATION_CHECK(A, B)
    if(!is_sor_convergent(A, w, step)) {
        PyErr_SetString(PyExc_ValueError, 
            "the SOR iteration is not convergent, "
            "set 'step' True to get more information");
        return NULL;
    }
    M = sor(A, B, w, delta, step);
    //M = sor_homework(A, B, w, delta, step);
    return (PyObject*)mym_matrix_create(M);
}

static PyObject* mym_pow_method(PyObject* self, PyObject* args, PyObject* kwds)
{
    Matrix *A;
    bool step = false;
    static char* kwlist[] = {"A", "step", NULL};
    if(!PyArg_ParseTupleAndKeywords(args, kwds, 
        "O&|O&:pow_method", kwlist, 
        &mym_matrix_converter, &A, 
        &mym_bool_converter, &step)) {
        return NULL;
    }
    return PyFloat_FromDouble(pow_method(A, step));
}

static PyObject* mym_norm(PyObject* self, PyObject* args, PyObject* kwds)
{
    Matrix *A;
    int p;
    bool step = false;
    static char* kwlist[] = {"A", "p", "step", NULL};
    if(!PyArg_ParseTupleAndKeywords(args, kwds, 
        "O&I|O&:norm", kwlist, 
        &mym_matrix_converter, &A,
        &p, &mym_bool_converter, &step)) {
        return NULL;
    }
    return PyFloat_FromDouble(norm(A, p, step));
}

//the only non-static function
PyMODINIT_FUNC PyInit_mym(void)
{
    PyObject* m;
    MYM_Matrix_Type.tp_new = PyType_GenericNew;

    if(PyType_Ready(&MYM_Matrix_Type) < 0)
        return NULL;
    m = PyModule_Create(&mym_module);
    if(NULL == m)
        return NULL;
    Py_INCREF(&MYM_Matrix_Type);
    PyModule_AddObject(m, "Matrix", (PyObject*)&MYM_Matrix_Type);
    return m;
}