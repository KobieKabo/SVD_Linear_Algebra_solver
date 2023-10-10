# SVD linear Algebra Solver
Contained in this repository is a rudimentary SVD solver for simple linear algebra problems. Additionally, is a spring-mass system script that uses the prior py script mentioned to solve for three different boundary condition problems.

# Running the Project

To run the scripts (detailed below):
```
$ python3 SVD_solver.py
OR
$ python3 Spring_mass.py
```
The terminal will then prompt you for input.


## Scripts
### `SVD_solver.py`
This script preforms SVD on a given matrix. Solving, A = U &Sigma; V<sup>T</sup> 

The SVD method decomposes a matrix into three other matrices: U, Sigma, and V. U is the left eigen-vector matrix, V is the right eigenvector matrix, and sigma is a diagonal matrix containing the square root of the eigen values of AA<sup>T</sup> and A<sup>T</sup>A, which should be equivalent. The script additionally prints a comparison to the numpy SVD blackbox solver.

There is one function call within this script and that is the following:

`svd_solver(A)`: Perform Singular Value Decomposition on matrix A; returning U, Sigma, and V, the condition numbers, and the inverse of A.

The script was written to take user input from the command line, including the size of the matrix A & the elements contained in it. However, the A inverse return has not been working and was something I was not able to figure out the issue for.\

**Example:**

```
PS C:\Users\Jakob\College\COE\COE 352\SVD_Linear_Algebra_solver> python3 .\SVD_Solver.py
Enter the number of rows in matrix A: 2
Enter the number of columns in matrix A: 2
Enter element A[1,1]: 1
Enter element A[1,2]: 0
Enter element A[2,1]: 0
Enter element A[2,2]: 1

A:
[[1. 0.]
 [0. 1.]]
___________________________

A Reconstruction:
[[1. 0.]
 [0. 1.]]
___________________________

U:
[[0. 1.]
 [1. 0.]]
___________________________

Sigma:
[[1. 0.]
 [0. 1.]]
___________________________

V:
[[0. 1.]
 [1. 0.]]
___________________________

Condition Numbers:
1.0
___________________________

A-Inverse:
[[0. 0.]
 [0. 0.]]
___________________________

Blackbox U:
[[1. 0.]
 [0. 1.]]
___________________________

BlackBox Sigma:
[1. 1.]
___________________________

BlackBox V:
[[1. 0.]
 [0. 1.]]
___________________________
```

I've written checks to catch when the inverse of a matrix does not exist; however, again like stated above the inverse calculation was not an issue I was able to resolve.

### `spring_mass.py`
The spring-mass script accepts user input for a spring mass system, asking for the boundary conditions, number of masses, spring constants and mass values, all by user input via the command line. The script then computes the equilibrium displacements **u**, the stiffness matrix **K**, in order to resovle the equatuon **f = Ku**, where **f** is the force vector.

```
PS C:\Users\Jakob\College\COE\COE 352\SVD_Linear_Algebra_solver> python3 .\Spring_mass.py
Please enter a boundary for the equation.
Enter 1 for one fixed end.
Enter 2 for two fixed ends.
Enter 3 for no fixed ends.

2
Please enter the number of masses

2
Enter the spring constant for spring  1 :
1
Enter the spring constant for spring  2 :
1
Enter the spring constant for spring  3 :
1
Enter the mass for mass  1 :
1
Enter the mass for mass  2 :
1
[[ 1.  0.]
 [-1.  1.]
 [ 0. -1.]]
[[1. 0. 0.]
 [0. 1. 0.]
 [0. 0. 1.]]
 
Your results are:
The mass displacements are:
 [[9.8]
 [9.8]]

The spring elongation is:
 [[ 9.80000000e+00]
 [ 3.55271368e-15]
 [-9.80000000e+00]]

The singular values and condition number for the A matrix are:
 [[1.73205081 0.        ]
 [0.         1.        ]
 [0.         0.        ]]        173234686.98031548
 
The singular values and condition number for the A matrix transposed are:
 [[1.73205081 0.         0.        ]
 [0.         1.         0.        ]]     1.7320508075688772

The singular values and condition number for the spring constant matrix are:
 [[1. 0. 0.]
 [0. 1. 0.]
 [0. 0. 1.]]     1.0

The K matrix is:
 [[ 2. -1.]
 [-1.  2.]]
 
The singular values and condition number for the K matrix are:
 [[3. 0.]
 [0. 1.]]        3.0
 ----------------------------------------
```
Note:
If the boundary condition of 3 is used, physically it makes no sense. As the masses would move without the springs, but for the solving of u shows that this is impossible. Thus the check for a two free end system returns an L2 value of infinity, indicating that our matrix is ill-conditioned.
