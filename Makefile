# Makefile for C++ project

# Compiler settings
CXX = g++
#CXXFLAGS = -march=native -std=c++20 -Wall -pthread -Ofast
#CXXFLAGS = -fopenmp -g -fsanitize=address -march=native -std=c++20 -w -pthread -Ofast
CXXFLAGS = -fopenmp -g -static-libasan -fsanitize=address -std=c++20 -w -O3 -march=native -I /home/vaidya2/boost_1_88_0
#CXXFLAGS = -march=native -std=c++20 -w -pthread -Ofast

# Source files
SRCS = annealer.cpp vertexing.cpp detanneal.cpp
HEADERS = annealer.hh vertexing.hh detanneal.hh

# Executable name
TARGET = annealer

# Default target
all: $(TARGET)

# Rule to create the executable
$(TARGET): $(SRCS) $(HEADERS)
	echo "Compiling $(SRCS)..."
	$(CXX) $(CXXFLAGS) -o $@ $(SRCS)
	echo "Build complete: $(TARGET)"

# Run the program 
run: $(TARGET)
	echo "Running $(TARGET)..."
	./$(TARGET) $(THREADS) $(STAGES) $(SAMPLES_PER_THREAD) $(SWEEPS) $(DA_SWEEPS) QPU_$(VERTICES)Vertices_$(TRACKS)Tracks_100Events/$(VERTICES)Vertices_$(TRACKS)Tracks_Event $(VERTICES)Vertices_$(TRACKS)Tracks /serializedEvents.json
	
	echo "Execution finished."

# Clean up generated files
clean:
	echo "Cleaning up..."
	rm -f $(TARGET)
	echo "Clean complete."

# Phony targets
.PHONY: all clean run
