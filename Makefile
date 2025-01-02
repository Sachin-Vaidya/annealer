# Makefile for C++ project

# Compiler settings
CXX = g++
CXXFLAGS = -std=c++17 -Wall -pthread -Ofast

# Source files
SRCS = annealer.cpp parser.cpp vertexing.cpp
HEADERS = parser.hh vertexing.hh ThreadPool.h

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
