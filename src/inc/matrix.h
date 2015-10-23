#ifndef _MATRIX_H_
#define _MATRIX_H_

#define MAX(x, y)   ((x) > (y) ? (x) : (y))
#define MIN(x, y)   ((x) < (y) ? (x) : (y))
#define DOUBLE_EQUAL(x, y)  ((x) - (y) < 1e-6 && (x) - (y) > -1e-6)
#define INF ((1LL << 31) - 1)
#define IS_SQUARE(M)    ((M)->m == (M)->n)
#define IS_VECTOR(V)    ((V)->n == 1)
#define IS_EQUATION(A, B)   (IS_SQUARE((A)) && IS_VECTOR((B)) && (A)->m == (B)->m)

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
static inline void swap_row(Matrix* M, int i, int j) {
    swap_dp(&M->mem[i], &M->mem[j]);
}

//创建n阶方阵
Matrix* create_matrix_n(int n);

//创建mxn阶方阵
Matrix* create_matrix(int m_, int n_);

//创建n阶单位阵
Matrix* create_eye(int n);

//销毁方阵
void destroy_matrix(Matrix* M);

void fill_matrix(Matrix* M, double x);

//矩阵相乘
Matrix* mul(Matrix* A, Matrix* B);

//计算矩阵的行列式
double det(Matrix* M);

//复制矩阵
Matrix* copy(Matrix* M);

//浅复制
Matrix* shallow_copy(Matrix* M);

//打印矩阵
void print_matrix(Matrix* M);
bool is_ordered_main_subdet(Matrix* M, Comp checker, double x);
bool is_symmetrical(Matrix* M);
bool is_const_similar(Matrix* M, Matrix* N, double* coeff);
Matrix* transpose(Matrix* M);
Matrix* inverse(Matrix* M);

void swap_matrix(Matrix* M, Matrix* N);

#endif