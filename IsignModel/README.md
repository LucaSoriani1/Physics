# Isign Model

<p>The Ising model is a physical model of a magnetic system consisting of a discrete $1D$ lattice of interacting spins that can be used to study ferromagnetism. In this model, there are $N$ quantum particles with spin $1/2$, where each spin can be oriented "upward" or "downward," represented by the values $+1$ and $-1$.

It has been used to study phase transitions in magnetic systems, where the system goes from an ordered state at low temperature to a disordered state at high temperature. Ising's model has also been used to describe a wide range of phenomena, including order in biological systems, protein aggregation, neural network formation, and even game theory. The energy of the system depends on the spin interactions and the energies of each particle, such that the Hamiltonian is:
$$
\widehat{H} = \lambda \sum_{j=1}^N \sigma^{(j)}_z + \sum_{j=1}^{N-1} \sigma^{(j)}_x \sigma^{(j+1)}_x,
$$
where $\sigma_x$ and $\sigma_z$ are two of the three Pauli matrices given by:
$$ \sigma_x = \begin{pmatrix}
0 & 1 \\
1 & 0
\end{pmatrix} \, \, \, \, \, \sigma_z = \begin{pmatrix}
 1 & 0 \\
0 & -1
\end{pmatrix} $$
The eigenvalues of the Ising model are distributed differently depending on the value of the parameter $\lambda$ representing the strength of the interaction between the spins.

In particular, for values of $\lambda$ close to 0, the eigenvalues become degenerate, that is, they cluster into a limited number of values, while for values of $\lambda$ greater than 0, the eigenvalues begin to distinguish themselves.

The energy of the system is given by the sum of the interaction energies between the spins. Specifically, each pair of neighboring spins contributes to the total energy of the system according to their interaction, which can be ferromagnetic or antiferromagnetic. In the former case, the spins tend to align parallel, while in the latter case they tend to align antiparallel.

Ising's model can be solved numerically using the Hamiltonian. The solution of the Hamiltonian makes it possible to calculate properties of the system, such as magnetization and free energy. Given such a Hamiltonian operator, the $j$-th particle contributes to that operator as:$$\widehat{H}_j = \lambda \left( \mathbb{I}_{2\times 2}^{(1)} \otimes \mathbb{I}_{2 \times 2}^{(2)} \otimes \ldots \otimes \sigma_z^{(j)} \otimes \ldots \otimes \mathbb{I}_{2\times 2}^{(N)}\right) \\ + \left( \mathbb{I}_{2 \times 2}^{(1)}  \otimes  \mathbb{I}_{2 \times 2}^{(2)} \otimes  \ldots \otimes \sigma_x^{(j)}\otimes \ldots \otimes \mathbb{I}_{2 \times 2}^{(N)}\right). $$ 

This script numerically solves the Ising model by calculating the Hamiltonian of the system. The Hamiltonian is represented as a $2^N \times 2^N$ matrix, where $N$ is the number of particles, and is calculated separately for the spin interactions along the $z$-axis and along the $x$-axis. Next, the eigenvalues and eigenvectors of the Hamiltonian matrix are calculated using the 'numpy.linalg.eigh' function. The code also includes a progress bar to monitor the progress of the calculation and an estimate of the expected time. 


The 'Matrix' class defines a square matrix initialized to $0$ of dimension 'dim'.

The 'pauli' function creates two Pauli matrices, 'sigma_x' and 'sigma_z', of size $2\times 2$.

The 'z_hamiltonian' function creates the z part of the Hamiltonian for the Ising model. The z Hamiltonian is a matrix of dimension $2^N \times 2^N$. The Pauli $z$ sigma_z matrix is used to create the $z$ Hamiltonian.

The 'x_hamiltonian' function creates the x part of the Hamiltonian repredented by a matrix of dimension $2^N \times 2N$. The Pauli $x$ sigma_x matrix is used to create the Hamiltonian $x$.

The 'expected_time' function calculates the expected execution time to complete a certain number of iterations based on the average time of previous iterations.
The 'progress_bar' function prints a progress bar with an estimate of the time remaining to complete the execution of the algorithm.

### Usage
The code can be run from the command line by executing the time_dependent_harmonic_oscillator() function. The user can adjust the following parameters:

* l_min: the minimum value of the parameter $\lambda$
* l_max: the maximum value of the parameter $\lambda$
* n_particles: the number of quantum particles in the system.

These parameters are set in their optimal configuration. You can change them and 'play' with them, but the computation may take a lot of time. As default will be plot the 2**n_particles eigevalues, but you can change it by setting the 'kappa' parameter

### Dependencies
This code requires the following packages:

* NumPy
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
