#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <stdlib.h>
#include <memory.h>

#define SIGN(x)    ((x) >= .0 ? 1 : -1)
#define MAX(x, y)   ((x) > (y) ? (x) : (y))
#define MIN(x, y)   ((x) < (y) ? (x) : (y))
#define DOUBLE_EQUAL(x, y)  ((x) - (y) < 1e-6 && (x) - (y) > -1e-6)
#define DOUBLE_EQUAL_DELTA(x, y, delta) ((x) - (y) < (delta) && (x) - (y) > -(delta))
#define INF ((1LL << 31) - 1)
#define IS_SQUARE(M)    ((M)->m == (M)->n)
#define IS_VECTOR(V)    ((V)->n == 1)
#define IS_EQUATION(A, B)   (IS_SQUARE((A)) && IS_VECTOR((B)) && (A)->m == (B)->m)
#define IS_SAME_SIZE(A, B) ((A)->m == (B)->m && (A)->n == (B)->n)
#define IS_VALID(M) (NULL != (M) && NULL != (M)->alloc && NULL != (M)->mem)


extern int step;

typedef enum{
    false = 0,
    true = 1
} bool;

typedef struct {
	int m;
	int n;
	double** mem;

    /*why we need this member 'alloc'?
    as member 'mem' is designed to make 
    line swap efficient. it possible that 
    swap(&mem[0], &mem[i]) can dirty the 
    memory allocation, and this is this 
    'alloc' for.*/
    double* alloc;
} Matrix;

typedef struct {
    int size;
    double* mem;
} Array;


typedef bool (*Comp)(double, double);
static inline bool gt(double x, double y) {
    return x > y;
}

static inline bool eq(double x, double y) {
    return DOUBLE_EQUAL(x, y);
}

static inline bool ne(double x, double y) {
    return !eq(x, y);
}

static inline bool lt(double x, double y) {
    return !gt(x, y) && ne(x, y);
}

void swap_dp(double** a, double** b);
void swap_d(double* a, double* b);
void swap_n(int* a, int* b);

static inline void swap_row(Matrix* M, int i, int j) {
    swap_dp(&M->mem[i], &M->mem[j]);
}

Array* create_array(int size);

void static inline destroy_array(Array** arr)
{
    if (NULL == arr || NULL == *arr)
        return;
    free((*arr)->mem);
    (*arr)->size = 0;
    *arr = NULL;
}

void static inline fill_array(Array* arr, double x)
{
    memset(arr->mem, x, sizeof(double) * arr->size);
}

//创建mxn阶方阵
Matrix* create_matrix(int m_, int n_);

/**
创建n阶方阵
*/
static inline Matrix* create_matrix_n(int n)
{
    return create_matrix(n, n);
}

static inline Matrix* create_vector(int m)
{
    return create_matrix(m, 1);
}

//创建n阶单位阵
Matrix* create_eye(int n);

Matrix* create_cons(int n, double x);

//销毁方阵
void destroy_matrix(Matrix** M);



void fill_matrix(Matrix* M, double x);

//矩阵相乘
Matrix* mul(const Matrix* A, const Matrix* B);
void mul_inp_L(Matrix* lhs, const Matrix* rhs);
void mul_inp_R(const Matrix* lhs, Matrix* rhs);
Matrix* mul_cons(const Matrix* A, double x);
void mul_cons_inp(Matrix* A, double x);
Matrix* sub(const Matrix* A, const Matrix* B);
void sub_inp(Matrix* A, const Matrix* B);
Matrix* add(const Matrix* A, const Matrix* B);
void add_inp(Matrix* A, const Matrix* B);
//计算矩阵的行列式
double det(Matrix* M);

//复制矩阵
Matrix* copy(const Matrix* M);
void copy_inp(Matrix* des, const Matrix* src);

void make_eye(Matrix* M);

//浅复制
Matrix* shallow_copy(Matrix* M);

//打印矩阵
void print_matrix(const Matrix* M);
bool is_ordered_main_subdet(Matrix* M, Comp checker, double x);
bool is_symmetrical(const Matrix* M);
bool is_const_similar(const Matrix* M, const Matrix* N, double delta, double* coeff);
Matrix* transpose(const Matrix* M);
void transpose_inp(Matrix* M);
Matrix* inverse(Matrix* M);

void swap_matrix(Matrix* M, Matrix* N);

unsigned long long factorial(size_t n);


#endif

