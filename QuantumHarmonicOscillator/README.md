# Quantum Harmonic Oscillator in 1D
The quantum harmonic oscillator is a mathematical model that describes the behavior of a particle oscillating around an equilibrium point constrained to move along an axis with a harmonic potential. The system is described by the Hamiltonian:
$$\hat{H} = \frac{\hat{p}^2}{2m} + \frac{1}{2}m\omega^2\hat{x}^2$$
where $m$ is the mass of the particle and $\omega$ is the frequency of the harmonic motion. 

The solution of the Schrödinger equation for the quantum harmonic oscillator is known and is given by the Hermite functions, which represent the stationary states of the oscillator. These states are quantized and their energy is given by the formula:
$$\psi_n(x) = \frac{1}{\sqrt{2^n n!}} \left(\frac{m\omega}{\pi \hbar}\right)^{\frac{1}{4}} e^{-\frac{m\omega x^2}{2 \hbar}} H_n \left(\sqrt{\frac{m\omega}{\hbar}}x \right). $$ 

where $n$ is the quantum number that identifies the energy state of the system, $\hbar$ is the reduced Planck constant, $H_n(x)$ are the Hermite polynomials of degree $n$. The eigenvalues of the system are the allowed energies of the particle in the harmonic potential, and the eigenfunctions are the wave functions describing the probability distribution of the particle in the harmonic potential.

The quantum harmonic oscillator has several applications in physics, particularly in quantum mechanics and quantum field theory. For example, it is used to describe the behavior of bound atoms and simple molecules, as well as to describe the behavior of subatomic particles such as electrons in a magnetic field. This is a very important mathematical model because it is one of the few quantum physical system models for which an exact analytical solution of the Schrödinger equation can be found. This makes the quantum harmonic oscillator one of the most studied models in quantum mechanics.
.

The code solves the eigenvalue problem of the quantum harmonic oscillator $1D$, setting $\hbar = \omega = 1$ for simplicity. In order to solve the problem, the discretization solution is adopted: in theory, the space is continuous, consisting of an infinite number of points. Unfortunately, the reasoning cannot be applied at the practical level where we can only deal with a finite number of points. So we want to adopt a very small step $dx$ in order to approximate the continuous model as best as possible. The discretization yields a simpler tridiagonal Hamiltonian $H_{discr}$, where in the pricipal diagonal we find the terms $-2+V(x_i)$ and in the sub-diagonals we find $-1$.

The 'quantum_harmonic_oscillatori function takes as input the minimum and maximum value of the $x$-axis, the discretization step of the $dx$ grid, and the number of eigenvalues to be computed.

First, the range of the $x$-axis, its mean value, and the number of steps within the range are calculated. Next, a matrix of the quantum harmonic oscillator potential and an associated diagonal matrix is created. The diagonal matrix consists of the sum of the potential and a constant. In addition, a subdiagonal to tridiagonal matrix with constant values is created, so that $H_{disc}$ is obtained. 

The eigenvalue problem is solved using the 'eigh_tridiagonal' function of the 'scipy.linalg' package. This function returns the eigenvalues and eigenvectors of the tridiagonal matrix from its three diagonals.

Finally, the 'print_results' function is used to print a table of the values of the expected eigenvalues, the calculated eigenvalues, and the percentage difference between the two. In addition, a graph is generated showing the eigenvalues and eigenvectors of the quantum harmonic oscillator.

### Usage
The code can be run from the command line by executing the main() function. The user can adjust the following parameters:

* x_min: the minimum value of the x-axis
* x_max: the maximum value of x
* discretization_dx: the step size of the discretization of the interval
* n_of_eigvs: number of eigenvalues to be calculated

These parameters are set in their optimal configuration. You can change them and 'play' with them, but the results may change and deviate from the theory. For example, changing the range (x_max - x_min) will allow you to calculate more correct eigenvalues 

The code outputs a table of the expected eigenvalues, the computed eigenvalues, and the percentage difference between the two. Additionally, a plot of the eigenvalues and eigenfunctions of the quantum harmonic oscillator is displayed.

### Dependencies
This code requires the following packages:

* NumPy
* SciPy
* Matplotlib
* tabulate

### Installation

This program requires Python 3 to be installed. All required packages are listed in the requirements.txt file. To install the dependencies, run:

```
pip install -r requirements.txt
```

### Fortran90 code

You will find the original code in Fortran90. To run the .f90 script, you will need a compiler, then some extension in Visual Studio code as:

- C/C++
- Modern Fortran
- Code Runner

Moreover, you need to install the 'lapack' library in order to run it.

### Modifications and Contributions

This program can be extended and modified to suit different purposes, such as modeling natural processes or generating random data for testing purposes. If you would like to contribute to the program or propose modifications, please send a pull request or contact the author.
