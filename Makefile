# Makefile for Neutron Diffusion Solver
CXX = g++
CXXFLAGS = -std=c++11 -Wall -O3
LDFLAGS =

# Directories
SRC_DIR = .
MATRIX_DIR = matrix_modules
LUP_DIR = lup
POINT_JACOBI_DIR = point_jacobi
GAUSS_SEIDEL_DIR = gauss_seidel
SOR_DIR = sor

# Include directories
INCLUDES = -I$(SRC_DIR) -I$(MATRIX_DIR) -I$(LUP_DIR) -I$(POINT_JACOBI_DIR) -I$(GAUSS_SEIDEL_DIR) -I$(SOR_DIR)

# Source files
SOURCES = $(SRC_DIR)/main.cpp \
		  $(MATRIX_DIR)/matrix_modules.cpp \
          $(LUP_DIR)/lup_solver.cpp \
          $(POINT_JACOBI_DIR)/point_jacobi_solver.cpp \
          $(GAUSS_SEIDEL_DIR)/gauss_seidel_solver.cpp \
          $(SOR_DIR)/sor_solver.cpp


# Object files
OBJECTS = $(SOURCES:.cpp=.o)

# Target executable
TARGET = diffusion_solver

# Default target
all: $(TARGET)

# Link object files to create executable
$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $(TARGET) $(LDFLAGS)

# Compile source files to object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Clean generated files
clean:
	rm -f $(OBJECTS) $(TARGET)
	rm -f output.txt

# Clean and rebuild
rebuild: clean all