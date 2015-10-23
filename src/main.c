#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <string.h>
#include "matrix.h"
#include "direct.h"
#include "norm.h"
#include "iterative.h"
#include "eigenvalue.h"

static char* app_name = "mym";

void usage(void) {
    printf("Usage:\n"
		   "    %s gauss <file> [step]      -use gauss elimination\n"
           "    %s megauss <file> [step]    -use main element gauss elimination\n"
           "    %s tri <file> [step]        -use triangular decomposition\n"
           "    %s metri <file> [step]      -use main element triangular decomposition\n"
           "    %s cholesky <file> [step]   -use cholesky decomposition\n"
           "    %s encholesky <file> [step] -use enhanced cholesky decomposition\n"
           "    %s charsing <file> [step]   -use charsing method\n"
           "    %s jacobi <file> [step]     -use jacobi iteration\n"
           "    %s gauss_seidel <file> [step] -use Gauss-Seidel iteration\n"
           "    %s sor <file> [step]        -use SOR iteration\n"
           "    %s help                     -show this usage\n"
           "\n", app_name, app_name, app_name, 
           app_name, app_name, app_name, app_name,
           app_name, app_name, app_name, app_name
        );
}

int main(int argc, char* argv[]) {
//int main(void) {
    //char* argv[3] = {"mym", "jacobi", "./inputs/m.txt"};
    //int argc = 3;
    int i, j, n;
    Matrix *A = NULL, *B = NULL, *X = NULL;
    if (argc < 2) {
	usage();	
	return -1;
    } else if(argc >= 2 && 0 == strcmp(argv[1], "help")) {
	usage();
	return -1;
    }
    
    app_name = argv[0];
    if(argc >= 3) {
		if (NULL == freopen(argv[2], "r", stdin)) {
			printf("error:can not find file:%s !\n", argv[2]);
			return -2;
		}
    } else {
	usage();
	return -1;
    }
    
    if (argc >= 4 && 0 == strcmp(argv[3], "step"))
        step = 1;

    scanf("%d", &n);
    A = create_matrix_n(n);
    B = create_matrix(n, 1);
    //X = create_matrix(n, 1);
    if(NULL == A || NULL == B) {
        printf("error: memory shortage!\n");
        return -1;
    }

    for(i = 0; i < A->m; i++)
        for(j = 0; j < A->n; j++)
            scanf("%lf", &A->mem[i][j]);
    for(i = 0; i < B->m; i++)
        scanf("%lf", &B->mem[i][0]);

	//printf("A det: %lf\n", det(A));
	printf("A matrix:\n");
	print_matrix(A);
	printf("\n");
	printf("B matrix:\n");
	print_matrix(B);
	printf("\n");

    if(0 == strcmp("gauss", argv[1])) {
        X = gauss_elim(A, B, step);
    } else if(0 == strcmp("megauss", argv[1])) {
        X = me_gauss_elim(A, B, step);
    } else if(0 == strcmp("tri", argv[1])) {
		X = tri_decomp(A, B, step);
    } else if(0 == strcmp("metri", argv[1])) {
		X = me_tri_decomp(A, B, step);
    } else if(0 == strcmp("cholesky", argv[1])) {
        X = cholesky_decomp(A, B, step);
    } else if( 0 == strcmp("encholesky", argv[1])) {
        X = en_cholesky_decomp(A, B, step);
    } else if(0 == strcmp("chasing", argv[1])) {
        X = chasing_method(A, B, step);
    } else if(0 == strcmp("jacobi", argv[1])) {
        X = jacobi_iter(A, B, 1e-6, step);
    } else if(0 == strcmp("gauss_seidel", argv[1])) {
        X = gauss_seidel_iter(A, B, 1e-6, step);
    } else if(0 == strcmp("sor", argv[1])) {
        X = sor(A, B, 1, 1e-6, step);
    } else if(0 == strcmp("pow", argv[1])) {
        double eigen_value = pow_method(A, step);
        printf("the biggest eigen value is: %lf\n", eigen_value);
    } else if(0 == strcmp("inv", argv[1])) {
        X = inverse(A);
    } else {
        /*
        #define N (5)
        Matrix* test = create_matrix(N, 1);
        const int arr[N] = {-10, 3, 5, -7, 8};
        int i;
        double no;
        for(i = 0; i < N; i++)
            test->mem[i][0] = arr[i];
        */
        Matrix* test = A;
        printf("norm: %lf\n", norm(test, 2, step));
    }

	if (NULL != X) {
		printf("the solution is:\n");
		print_matrix(X);
	}
	destroy_matrix(A);
	destroy_matrix(B);
	destroy_matrix(X);

    return 0;
}