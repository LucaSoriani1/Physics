# Quantum Harmonic Oscillator in 1D
The code solves the eigenvalue problem of the 1D quantum harmonic oscillator, which is a theoretical model of a physical system consisting of a particle constrained to move along an axis with a harmonic potential. The system is described by the Hamiltonian:

$$\hat{H} = \frac{\hat{p}^2}{2m} + \frac{1}{2}m\omega^2\hat{x}^2$$

where $\hat{p}$ and $\hat{x}$ are the momentum and position operators, respectively, $m$ is the mass of the particle, and $\omega$ is the frequency of the harmonic motion.

The eigenvalues of the system are the allowed energies of the particle in the harmonic potential, while the eigenfunctions are the wave functions that describe the probability distribution of the particle in the harmonic potential.

The eigenvalues and eigenfunctions of the system are given by:

$$ E_n = \hbar \omega \left( n + \frac{1}{2} \right) $$

$$ \psi_n(x) = \frac{1}{\sqrt{2^n n!}} \left(\frac{m\omega}{\pi \hbar}\right)^{\frac{1}{4}} e^{-\frac{m\omega x^2}{2 \hbar}} H_n \left(\sqrt{\frac{m\omega}{\hbar}}x\right) $$

where $n$ is the quantum number that identifies the energy state of the system, $\hbar$ is the reduced Planck constant, $H_n(x)$ are the Hermite polynomials of degree $n$, and $\omega$ is the angular frequency of the harmonic motion.

In the code, we set $\hbar = \omega =  1$ and the solution of the eigenvalues is obtained through the numerical solution of the time-dependent Schrödinger equation, which is rewritten as a position-dependent eigenvalue problem. In particular, a reduced-banded symmetric tridiagonal matrix representing the harmonic potential is created, and the eigenvalues and eigenfunctions are calculated using the eigh_tridiagonal function of the scipy.linalg module. Finally, the graph of the eigenfunctions and eigenvalues of the system is displayed.

In summary, the code implements the numerical solution of the eigenvalue problem for the 1D quantum harmonic oscillator. The quantum harmonic oscillator is a quantum mechanical system that describes the behavior of a particle constrained by a harmonic potential, i.e., a particle subject to a force proportional to its distance from the equilibrium position. In particular, the quantum harmonic oscillator is a very important model in quantum mechanics because it provides an example of a many-body system that can be exactly solved and it also appears in many other physical contexts.

The code implements the numerical algorithm to calculate the energies and wave functions of the 1D quantum harmonic oscillator. In particular, the code defines the system interval, step size, and the number of eigenvalues and eigenvectors to calculate. Subsequently, the code defines the potential energy matrix of the quantum harmonic oscillator and uses the eigh_tridiagonal function of the scipy.linalg module to solve the eigenvalue problem for the tridiagonal matrix associated with the potential energy matrix. Finally, the code displays the energies and wave functions of the 1D quantum harmonic oscillator.

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
