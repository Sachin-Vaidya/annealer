# Makefile for C++ project

# Compiler settings
CXX = g++
CXXFLAGS = -march=native -std=c++17 -Wall -pthread -Ofast

# Source files
SRCS = annealer.cpp vertexing.cpp detanneal.cpp
HEADERS = annealer.hh vertexing.hh detanneal.hh

# Executable name
TARGET = annealer

# Default target
all: $(TARGET)

# Rule to create the executable
$(TARGET): $(SRCS) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ $(SRCS)

# Run the program
run: $(TARGET)
	./$(TARGET)

# Clean up generated files
clean:
	rm -f $(TARGET)

# Phony targets
.PHONY: all clean run
