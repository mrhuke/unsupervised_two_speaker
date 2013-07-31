# Variable definition
CC = g++
LD = g++
CFLAGS = 
LDFLAGS = -lm 
PROG = segmentation
SRCS = filter.cpp tool.cpp gammaTone.cpp segmentation.cpp sgementation.cpp mScaleInten.cpp 
OBJS = filter.o tool.o gammaTone.o segmentation.o sgementation.o mScaleInten.o
RM = rm
# end

all:$(PROG)

# Link
 $(PROG):$(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) -o $(PROG)

# Compile
%.o:%.c
	$(CC) $(CFLAGS) -c $<

clean:
	$(RM) $(OBJS)

