# Makefile for Neutron Diffusion Solver
CXX = g++
MPICXX = mpicxx
CXXFLAGS = -std=c++11 -Wall -O3
LDFLAGS =

# Directories
SRC_DIR = .
MATRIX_DIR = matrix_modules
LUP_DIR = lup
POINT_JACOBI_DIR = point_jacobi
GAUSS_SEIDEL_DIR = gauss_seidel
SOR_DIR = sor
PARALLEL_PJ_DIR = parallel_pj
PARALLEL_GS_DIR = parallel_gs

# Include directories
INCLUDES = -I$(SRC_DIR) -I$(MATRIX_DIR) -I$(LUP_DIR) -I$(POINT_JACOBI_DIR) \
           -I$(GAUSS_SEIDEL_DIR) -I$(SOR_DIR) -I$(PARALLEL_PJ_DIR) -I$(PARALLEL_GS_DIR)

# Source files
SOURCES = $(SRC_DIR)/main.cpp \
          $(MATRIX_DIR)/matrix_modules.cpp \
          $(LUP_DIR)/lup_solver.cpp \
          $(POINT_JACOBI_DIR)/point_jacobi_solver.cpp \
          $(GAUSS_SEIDEL_DIR)/gauss_seidel_solver.cpp \
          $(SOR_DIR)/sor_solver.cpp \
          $(PARALLEL_PJ_DIR)/parallel_pj_solver.cpp \
          $(PARALLEL_GS_DIR)/parallel_gs_solver.cpp

# Object files
OBJECTS = $(SOURCES:.cpp=.o)

# Target executable
TARGET = diffusion_solver

# Default target (serial version)
all: serial

# Serial version
serial: COMPILER = $(CXX)
serial: $(TARGET)

# Parallel version with MPI
parallel: COMPILER = $(MPICXX)
parallel: CXXFLAGS += -DUSE_MPI
parallel: $(TARGET)

# Link object files to create executable
$(TARGET): $(OBJECTS)
	$(COMPILER) $(CXXFLAGS) $(OBJECTS) -o $(TARGET) $(LDFLAGS)

# Compile source files to object files
%.o: %.cpp
	$(COMPILER) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Clean generated files
clean:
	rm -f $(OBJECTS) $(TARGET)
	rm -f output.txt
	rm -f Flux

# Clean and rebuild serial version
rebuild: clean all

# Clean and rebuild parallel version
rebuild-parallel: clean parallel

.PHONY: all serial parallel clean rebuild rebuild-parallel