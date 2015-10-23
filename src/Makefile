CFLAGS = -Wall
SRCS = main.c matrix.c direct.c eigenvalue.c norm.c iterative.c

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

