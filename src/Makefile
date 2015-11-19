CFLAGS = -Wall
SRCS = main.c matrix.c poly.c direct.c nonlinear.c approximation.c \
	eigenvalue.c norm.c iterative.c interpolation.c integration.c

TARGET = mym

OBJS = $(SRCS:.c=.o)

$(TARGET) : $(OBJS)
	cc -o $@ $^ -lm
INC = -I./inc
%.o : %.c
	cc $(CFLAGS) $(INC) -o $@ -c $<
	
.PHONY : clean
clean :
	-rm -rf $(TARGET) $(OBJS)

