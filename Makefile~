objs = main.o matrix.o
target = mym
$(target) : $(objs)
	cc -o $(target) $(objs)
main.o : main.h
matrix.o : matrix.h

.PHONY : clean
clean :
	-rm $(target) $(objs)

