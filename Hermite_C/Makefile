# Compiler
CC = g++

# Compiler flags
CFLAGS = -std=c++11 -Wall 

# Linker flags
LDFLAGS = -larmadillo

# Target executable
TARGET = hermite_fem

# Source files
SRCS = hermite_fem.cpp
# Header files
HDRS = hermite_bases.h hermite_polynomials.h

# Object files
OBJS = $(SRCS:.cpp=.o)

# Rule to compile .cpp files to .o files
%.o: %.cpp $(HDRS)
	$(CC) $(CFLAGS) -c $< -o $@

# Rule to build the executable
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $(TARGET) $(LDFLAGS)

# Rule to clean the object files and executable
clean:
	rm -f $(OBJS) $(TARGET)
