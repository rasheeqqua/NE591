# Inlab8: Neutron Transport Equation Solver

This program implements a solver for the discrete ordinates approximation of the steady state, one-speed neutron
transport equation in one-dimensional slab geometry. It uses the source iteration algorithm to solve for the
neutron flux in a homogeneous medium with vacuum boundary conditions.

## Compilation

To compile the code, use the following command:

```bash
chmod +x compile.sh
./compile.sh
```

## Execution

To run the program, use the following command:

```bash
bsub < submit.sh
```

To read the output file, use the following command:
```bash
cat output.txt
```

Where:
- `input.txt` is the input file containing problem parameters
- `output.txt` is the file where results will be written

## Status

Status: Operational

## Source Code Files

- `main.cpp` - Contains the full implementation of the neutron transport solver including the source iteration algorithm

## Problem Description

### Physical Problem

This program solves the steady state, one-speed neutron transport equation in one-dimensional slab geometry.
The problem models neutron transport through a homogeneous medium of thickness L with vacuum boundary conditions
on both edges of the slab.

The continuum equation being solved is:

```
μ(∂ψ(x,μ)/∂x) + σₜψ(x,μ) = σₛΦ(x) + q
```

with boundary conditions:
- ψ(0,μ) = 0 for μ > 0
- ψ(L,μ) = 0 for μ < 0

Where:
- μ is the cosine of the polar angle
- ψ(x,μ) is the angular neutron flux
- Φ(x) is the scalar neutron flux
- σₜ is the total macroscopic cross section
- σₛ is the scattering macroscopic cross section
- q is a uniformly distributed, isotropic fixed source of neutrons

### Numerical Method

The program uses the discrete ordinates (SN) approximation with the Diamond Difference method for spatial discretization.
The source iteration algorithm is implemented to handle the scattering source term:

1. Initialize the scalar flux to zero
2. Compute the scattering source using the previous iteration's scalar flux
3. Sweep the grid along positive μ directions (left to right)
4. Sweep the grid along negative μ directions (right to left)
5. Update the scalar flux and check for convergence
6. Repeat steps 2-5 until convergence or maximum iterations reached

### Input

The input file should contain the following parameters in order:
- N: Number of angles in the quadrature set (the total number of angles will be 2N)
- I: Number of spatial cells
- σₜ: Total macroscopic cross section
- σₛ: Scattering macroscopic cross section
- q: Fixed source strength
- L: Slab width
- ε: Convergence tolerance
- kₘₐₓ: Maximum number of iterations
- Quadrature weights and points (ωₙ μₙ pairs, one per line for n=1,...,N)

Example input file:
```
1 10
1.0 0.5
1.0 1.0
1E-5 100
0.5 0.57735
```

### Output

The program outputs:
- Problem description and input parameters
- Convergence status and number of iterations
- Final relative error
- Scalar flux values for each spatial cell
- Execution time

### Limitations

- The program assumes a homogeneous medium (constant cross sections)
- Only vacuum boundary conditions are implemented
- The source is assumed to be uniform throughout the domain
- Limited to 1D slab geometry