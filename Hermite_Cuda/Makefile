# Compiler
CC = nvcc

# Compiler flags
CFLAGS = -std=c++11 -Xcompiler -Wall 

# Linker flags
LDFLAGS = -larmadillo

# Target executable
TARGET = hermite_fem

# Source files
SRCS = hermite_fem.cuh
# CUDA source files
CU_SRCS = hermite_bases.cu
# Header files
HDRS = hermite_bases.cuh hermite_polynomials.cuh

# Object files
OBJS = $(SRCS:.cpp=.o)
# CUDA object files
CU_OBJS = $(CU_SRCS:.cu=.o)

# Rule to compile .cpp files to .o files
%.o: %.cpp $(HDRS)
	$(CC) $(CFLAGS) -c $< -o $@

# Rule to compile .cu files to .o files
%.o: %.cu $(HDRS)
	$(CC) $(CFLAGS) -c $< -o $@

# Rule to build the executable
$(TARGET): $(OBJS) $(CU_OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(CU_OBJS) -o $(TARGET) $(LDFLAGS)

# Rule to clean the object files and executable
clean:
	rm -f $(OBJS) $(CU_OBJS) $(TARGET)
