# Makefile for C++ project

# Compiler settings
CXX = g++
#CXXFLAGS = -march=native -std=c++14 -Wall -pthread -Ofast
#CXXFLAGS = -fopenmp -g -fsanitize=address -march=native -std=c++14 -w -pthread -Ofast
CXXFLAGS = -fopenmp -g -std=c++14 -w -O2 -march=native
#CXXFLAGS = -march=native -std=c++14 -w -pthread -Ofast

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
	#./$(TARGET) $(THREADS) $(STAGES) $(SAMPLES_PER_THREAD) $(SWEEPS) 5V_30T/events_5V_30T_ 5Vertices_30TracksPerVertex .json
	./$(TARGET) $(THREADS) $(STAGES) $(SAMPLES_PER_THREAD) $(SWEEPS) QPU_3Vertices_15TracksTotal_100Events/3Vertices_15Tracks_Event 3Vertices_5TracksPerVertex /serializedEvents.json
	echo "Execution finished."

# Clean up generated files
clean:
	echo "Cleaning up..."
	rm -f $(TARGET)
	echo "Clean complete."

# Phony targets
.PHONY: all clean run
