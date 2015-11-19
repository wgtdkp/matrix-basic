#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <string.h>
#include "matrix.h"
#include "direct.h"
#include "norm.h"
#include "iterative.h"
#include "eigenvalue.h"
#include "poly.h"
#include "interpolation.h"
#include "nonlinear.h"
#include "integration.h"
#include "approximation.h"


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

static double foo(double x, void* coeff)
{
    return 1.0 / x / x;
}

static double bar(double x, void* coeff)
{
    return x * x * exp(-x);
}


int main(int argc, char* argv[]) {
//int main(void) {
    //char* argv[3] = {"mym", "qr", "./inputs/m.txt"};
    //int argc = 3;
    int i, j, n;
    Matrix *A = NULL, *B = NULL, *X = NULL;
    if (argc <= 2) {
    	usage();	
    	return -1;
    }
    
    app_name = argv[0];
	if (NULL == freopen(argv[2], "r", stdin)) {
		printf("error:can not find file:%s !\n", argv[2]);
		return -2;
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
        double eigen_value = pow_method(A, 1e-6, step);
        printf("the biggest eigen value is: %lf\n", eigen_value);
    } else if(0 == strcmp("inv", argv[1])) {
        X = inverse(A);
    } else if (0 == strcmp("qr", argv[1])) {
        //Matrix* Q = qr_decomp_G(A, true);
        //Matrix* R = A;
        //printf("Q matrix: \n");
        //print_matrix(Q);
        //to_hessenberg(A, true);   
        int i;
        Complex* eigenvals;
        eigenvals = (Complex*)malloc(sizeof(Complex) * 2 * (A->m - 1));
        qr_method(A, 1e-3, true);
        eigenvalues(eigenvals, A);
        printf("eigenvalues: \n");
        for (i = 0; i < A->m; i++) {
            print_complex(&eigenvals[i]);
            printf("\n");
        }
    } else if (0 == strcmp("poly", argv[1])) {
        Poly* poly = NULL;
        Poly* tmp = NULL;
        //poly = create_poly(3, 0);
        //printf("%.8lf\n", poly_value(poly, 2));
        //poly_add_inp(&poly, create_poly(1, 3));
        //poly_add_inp(&poly, create_poly(-3, 2));
        //poly_add_inp(&poly, create_poly(2, 1));
        //print_poly(poly);
        //poly_add_inp(&tmp, create_poly(1, 1));
        //poly_add_inp(&tmp, create_poly(-1, 0));
        //tmp = poly_copy(poly);
        //tmp = create_poly(1, 0);
        //poly_pow_inp(&poly, 3);
        //print_poly(poly);
        //printf("poly_value: %.6lf\n", poly_value(poly, 5));
        //poly_diff_inp(&poly);
        //poly_div_inp(&poly, tmp);
        //print_poly(poly);
        //destroy_poly(&poly);
        //poly_div_inp(&poly, tmp);
        //print_poly(poly);
    } else if (0 == strcmp("ip", argv[1])) {
        int i;
        Poly* poly;
        Poly** sp;
        #define n (6)
        double xarr[n] = {0, 1, 2, 3, 4, 5};
        double yarr[n] = {1, 2, 1, 2, 1, 3};
        //poly = newton(xarr, yarr, n);
        //poly = lagrange(xarr, yarr, n);
        sp = spline(xarr, yarr, n, 0, 0, 0);
        for (i = 0; i < n - 1; i ++) {
            printf("between [%.4lf, %.4lf]: ");
            print_poly(sp[i]);
            printf("diff at %.4lf is %.8lf\n", xarr[i], poly_value(sp[i], xarr[i]) - yarr[i]);
            printf("diff at %.4lf is %.8lf\n", xarr[i+1], poly_value(sp[i], xarr[i+1]) - yarr[i+1]);
        }
        //print_poly(poly);
        //for (i = 0; i < n; i++) {
        //    double value = poly_value(poly, xarr[i]);
        //    printf("poly_value(%.8lf) - yarr[%d] = %.8lf\n", xarr[i], i, value - yarr[i]);
        //}
    } else if (0 == strcmp("quadratic", argv[1])) {
        Complex* roots;
        Poly* poly = NULL;
        poly_add_inp(&poly, create_poly(1, 2));
        poly_add_inp(&poly, create_poly(6.94568726, 1));
        poly_add_inp(&poly, create_poly(-1222.34529948, 0));
        roots = (Complex*)malloc(sizeof(Complex) * 2 * (A->m - 1));
        quadratic(roots, poly);
        print_poly(poly);
        print_complex(&roots[0]);
        printf("\n");
        print_complex(&roots[1]);
        printf("\n");
    } else if (0 == strcmp("romberg", argv[1])) {
        double res;
        res = romberg(foo, 0.2, 1, 1e-5);
        printf("integration: %lf\n", res);
    } else if (0 == strcmp("simpson", argv[1])) {
        double res;
        res = simpson(foo, 0.2, 1, 1024);
        printf("integration: %lf\n", res);
    } else if (0 == strcmp("integration", argv[1])) {
        double res;
        res = integration(foo, 0.2, 1, 0.00133333);
        printf("integration: %lf\n", res);
    } else if (0 == strcmp("diff", argv[1])) {
        double res;
        //printf("%lf\n", exp(1));
        res = diff(bar, 0.5, 1e-12);
        printf("diff: %.13lf\n", res);
    } else {
        Matrix* test = A;
        printf("norm: %lf\n", norm(test, 2, step));
    }

	if (NULL != X) {
		printf("the solution is:\n");
		print_matrix(X);
	}
	destroy_matrix(&A);
	destroy_matrix(&B);
	destroy_matrix(&X);

    return 0;
}
