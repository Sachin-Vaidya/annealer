# Makefile for C++ project

# Compiler settings
CXX = g++
#CXXFLAGS = -march=native -std=c++20 -Wall -pthread -Ofast
#CXXFLAGS = -fopenmp -g -fsanitize=address -march=native -std=c++20 -w -pthread -Ofast
CXXFLAGS = -fopenmp -g -fsanitize=address -std=c++20 -w -O3 -march=native -I /home/vaidya2/Bell/boost_1_88_0
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
	#./$(TARGET) $(THREADS) $(STAGES) $(SAMPLES_PER_THREAD) $(SWEEPS) $(DA_SWEEPS) 5V_30T/events_5V_30T_ 5Vertices_30TracksPerVertex .json
	#./$(TARGET) $(THREADS) $(STAGES) $(SAMPLES_PER_THREAD) $(SWEEPS) $(DA_SWEEPS) QPU_3Vertices_15Tracks_100Events/3Vertices_15Tracks_Event 3Vertices_15Tracks /serializedEvents.json
	./$(TARGET) $(THREADS) $(STAGES) $(SAMPLES_PER_THREAD) $(SWEEPS) $(DA_SWEEPS) QPU_4Vertices_20Tracks_100Events/4Vertices_20Tracks_Event 4Vertices_20Tracks /serializedEvents.json
	#./$(TARGET) $(THREADS) $(STAGES) $(SAMPLES_PER_THREAD) $(SWEEPS) $(DA_SWEEPS) QPU_2Vertices_10Tracks_100Events/2Vertices_10Tracks_Event 2Vertices_10Tracks /serializedEvents.json
	
	echo "Execution finished."

# Clean up generated files
clean:
	echo "Cleaning up..."
	rm -f $(TARGET)
	echo "Clean complete."

# Phony targets
.PHONY: all clean run
