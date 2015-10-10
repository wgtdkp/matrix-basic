objs = main.o matrix.o
target = mym
$(target) : $(objs)
	cc -o $(target) $(objs) -lm
main.o : main.h
matrix.o : matrix.c matrix.h
	cc -c  matrix.c -lm

.PHONY : clean
clean :
	-rm $(target) $(objs)

