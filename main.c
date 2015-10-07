#include "matrix.h"
#include <stdio.h>
#include <memory.h>
#include <string.h>

static char* app_name = "mym";

void usage(void) {
    printf("Usage:\n"
		   "    %s gauss file [step]    -use gauss elimination\n"
           "    %s megauss file [step]  -use main element gauss elimination\n"
           "    %s tri file [step]      -use triangular decomposition\n"
           "    %s metri file [step]    -use main element triangular decomposition\n"
           "\n"
        );
}

int main(int argc, char* argv[]) {
//int main(void) {
    //char* argv[3] = {"mym", "metri", "stdin.txt"};
    //int argc = 3;
    int i, j, n;
    Matrix *A = NULL, *B = NULL, *X = NULL;
	if (argc < 2)
		return -1;

    if(argc >= 3)
		if (NULL == freopen(argv[2], "r", stdin)) {
			printf("error:can not find file:%s !\n", argv[2]);
			return -2;
		}
    
	if (argc >= 4 && 0 == strcmp(argv[3], "step"))
		print_steps = 1;

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
        X = gauss_elim(A, B);
    } else if(0 == strcmp("megauss", argv[1])) {
        X = me_gauss_elim(A, B);
    } else if(0 == strcmp("tri", argv[1])) {
		X = tri_decomp(A, B);
    } else if(0 == strcmp("metri", argv[1])) {
		X = me_tri_decomp(A, B);
    }

	if (NULL != X) {
		printf("the solution is:\n");
		print_matrix(X);
	}
	free(A);
	free(B);
	free(X);

    return 0;
}