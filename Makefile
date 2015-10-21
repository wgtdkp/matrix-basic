objs = main.o matrix.o direct.o norm.o iterative.o eigenvalue.o

target = mym

$(target) : $(objs)
	cc -o $(target) $(objs) -lm

main.o : main.c main.h
	cc -c  main.c -Wall

matrix.o : matrix.c matrix.h
	cc -c  matrix.c -Wall

direct.o : direct.c direct.h
	cc -c direct.c -Wall

norm.o : norm.c norm.h
	cc -c norm.c -Wall

iterative.o : iterative.c iterative.h
	cc -c iterative.c -Wall

eigenvalue.o : eigenvalue.c eigenvalue.h
	cc -c eigenvalue.c -Wall


.PHONY : clean
clean :
	-rm $(target) $(objs)

