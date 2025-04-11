### Nonlinear Neutron Diffusion Equation Solver

#### Compilation
Go in the source directory containing all the files and type this in terminal:
```bash
make
```

#### Execution
```bash
cd scripts
bsub < submit.sh
```

#### Status
- Operational

#### Source Code
- The mathematical implementation of the fixed-point iteration algorithm is located at:
    `./fpi/fixed_point_iteration.cpp`

#### Problem Description
This program solves a nonlinear version of the neutron diffusion equation using Fixed-Point Iterations.
The nonlinearity comes from a temperature-dependent removal term where ρ = ρ₀ + β/√φ.

#### Input
Input file "input.txt" containing:
- Stopping criterion and maximum iterations 
- Number of nodes per dimension 
- ρ₀ and β (linear and nonlinear removal terms)
- L (domain side-length)



Output

"output.txt": Contains calculation summary, including convergence information
"Flux": Contains the final flux values for each node in format "i j flux"

Limitations

Assumes vacuum boundary conditions
Uses a simple first-order fixed-point iteration method
Convergence may be slow for certain parameter combinations
Very small flux values may affect convergence due to the 1/√φ term