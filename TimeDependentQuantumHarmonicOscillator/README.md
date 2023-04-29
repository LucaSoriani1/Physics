# Time Dependent Quantum Harmonic Oscillator
The system is described by a Hamiltonian that depends on time as shown in equation $$
\widehat{H} = \frac{\widehat{p}^{\,2}}{2} + \frac{\left(\widehat{q} - q_0(t)\right)^2}{2},$$ where $\widehat{p}$ and $\widehat{q}$ are the momentum and position operators respectively, and $q_0(t) = t/T$ where $t \in [0, T]$.

The ground state $\ket{\psi_0}$ of the time-independent harmonic oscillator is evolved in time using the split order method. For a time-independent Hamiltonian, the ground state $\ket{\psi(t)}$ can be obtained by applying the time evolution operator $e^{-i\widehat{H}t}$ to the initial state $\ket{\psi_0}$, such that:
$$
\ket{\psi(t)} = e^{-i\widehat{H}t} \ket{\psi_0}.
$$
I can
use the fact that, if I take a small enough time-step $dt$ in order to solve the system: I discretize the system. I can evolve the state using
$$
\ket{\psi(t+dt)} \approx e^{-i\widehat{H}\dd t} \ket{\psi(t)}$$

 However, for a time-dependent Hamiltonian, the time evolution operator cannot be computed directly, and the split order method is used instead.

The split order method involves approximating the evolution operator $e^{-i\widehat{H} dt}$ using the Baker-Campbell-Hausdorff formula. This approximation can be written as $$
e^{-i\widehat{H} dt} \simeq e^{-i\frac{\widehat{V}}{2} dt} e^{-i\frac{\widehat{T}}{2} dt} e^{-i\frac{\widehat{V}}{2} dt},$$ where $\widehat{T}=\widehat{p}^{,2}/2$ is the kinetic energy operator and $\widehat{V}=(\widehat{q}-q_0(t))^2/2$ is the potential energy operator.

The time evolution of the ground state can be computed using equation 
$$
\psi(x,t+ dt) \simeq e^{-i\frac{\widehat{V}}{2} dt} \mathcal{F}^{-1} \left[ e^{-i\frac{p^2}{2} dt} \mathcal{F}\left[ e^{-i\frac{\widehat{V}}{2}  dt} \psi(x,t)\right] \right],
$$
 which involves a Fourier transform to obtain an eigenstate in the momentum representation, followed by an anti-Fourier transform to return to the position representation. The time evolution is performed iteratively by taking small time steps $ dt$ up to a maximum time $T$. 

To solve the Schrödinger equation for this Hamiltonian, one can discretize the $x$ space and write the Hamiltonian as a tridiagonal matrix. The discretization of space is performed by dividing the spatial interval $[-L, L]$ into $N$ equal parts of width $dx$, where $L$ is the size of the system and $N$ is the number of grid points.

The discretization of $x$ space implies that the position $x_i$ is given by $x_i = i*dx$ where $i = 0, 1, 2, ..., N$.

In this discretization, eache component of the Hamiltonian can be approximated by the following tridiagonal matrix:

$$
H_{ij} = \left[\left(-\frac{m}{2 dx^2}\right) + \frac{m}{2} w^2 x_i^2\right] \delta_{ij} - \frac{m}{2 dx^2}  \delta_{i,j-1} - \frac{m}{2 dx^2} \delta_{i,j+1}
$$

where $\delta_{ij}$ is the Kronecker delta, which is equal to $1$ if $i = j$ and $0$ otherwise.

This code simulates the time evolution of a quantum harmonic oscillator. The simulation is performed by numerically solving the time-dependent Schrödinger equation using the split order method and discretization. 

The code uses the 'NumPy' library to work with matrices and vectors, the 'SciPy' library to compute the eigenvalues and eigenvectors of a tridiagonal matrix, and the 'Matplotlib' library to create graphs.

The 'tridiagonal_and_eigenproblem' function creates the tridiagonal Hamiltonian matrix and calculates its eigenvalues and eigenvectors using SciPy's 'eigh_tridiagonal' function.

The 'momentum_space' function calculates the moment space representation of the wave function obtained from the eigenvector of the Hamiltonian. The function normalizes the wave function and calculates the values of the k moments, which are needed for the momentum space representation. 

The 'plot' function generates a graph of the wave function in space and time. Specifically, the function plots the wave function in space, the wave function in moment space, and the probability density in space and time.

The code includes several functions for calculating the expected time to complete a given number of iterations, displaying a progress bar with an estimate of the time remaining, normalizing a vector, and plotting the evolution of the wave function over time and the average position of the particle.

### Usage
The code can be run from the command line by executing the time_dependent_harmonic_oscillator() function. The user can adjust the following parameters:

* min: the minimum value of the x-axis
* max: the maximum value of x
* dx: the step size of the discretization for the space interval
* dt: the step size of the discretization for the time interval
* dk: the step size of the discretization for the momentum space
* quantum_number: the # of the eigenvalue/eigenvector computed
* T: total time interval

These parameters are set in their optimal configuration. You can change them and 'play' with them, but the results may change and deviate from the theory or the computation may take a lot of time. As default the first eigenvalue will be computed you can change it by varying the 'quantum_number'.  

The code outputs will be two images: the first will be the plot in the time of the eigenfunction vs the potential; the second one will be the mean position of the particle in the potential lie and the minimum of the potential and their variations in time.

### Dependencies
This code requires the following packages:

* NumPy
* SciPy
* Matplotlib
* time

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

Moreover, you need to install the 'lapack' library and 'Gnuplot'(adding it to the PATH variable) in order to run it.

### Modifications and Contributions

This program can be extended and modified to suit different purposes, such as modeling natural processes or generating random data for testing purposes. If you would like to contribute to the program or propose modifications, please send a pull request or contact the author.
