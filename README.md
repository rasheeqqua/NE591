# In-Lab Assignment 2

#### a) To compile the source code, go into the decompressed `inlab2` folder and open a terminal there. Then type in the following terminal commands:
```mkdir build```
```cd build```
```cmake -S .. -B .```
```make```

#### b) To run the code type this command in the terminal:
```./inlab2```

#### c) The code is: `operational`.

#### d) Brief Description of the Solved Problem
The program reads => n: number of interpolation points, m: number of evaluation points.
The user can either: Provide a file containing (x,y) pairs (Option 1)
Input x-values and compute y by calling the function f(x) (Option 2).
After reading the input data, the code prints the interpolation interval [a,b].
It then creates m equally spaced evaluation points in [a,b].
Finally, it prints a table of results, including placeholders for interpolation and error (which will be implemented in Outlab 2).

Variable declaration document:
n => Number of interpolation points
m => Number of evaluation points
userX => Stores the n user-provided x-values for interpolation points
choice => Chooses Option 1 (read from file) or Option 2 (evaluate using f(x))
xy => Stores(x,y) pairs (either read from file or computed by evaluate)
a, b => Defines the interpolation interval [a,b]
xEval => Stores the m equally spaced evaluation points in [a,b]

Algorithm Steps: 
Step #1: Print introductory header (assignment name, author, date).
Step #2: Read integers n and m.
Step #3: Validate n and m.
Step #4: Read user-provided interpolation points userX.
Step #5: Prompt and read the user’s choice (Option 1: file, Option 2: function).
Step #6: If choice == 1, read (x, y) pairs from a file into xy. If choice == 2, compute y using evaluate() for each userX[i] and populate xy.
Step #7: Print a confirmation of the input, plus the interval [a,b].
Step #8: Generate m evaluation points in [a,b] by calling createEvaluationPoints().
Step #9: Print the final table of results, with placeholders for interpolation and error.
Step #10: End the program with a completion message.
