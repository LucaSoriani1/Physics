# Isign Model

The Ising model is a physical model of a magnetic system consisting of a 1D discrete lattice of interacting spins that can be used to study ferromagnetism. In this model, there are N quantum particles with spin $1/2$, where each spin can be oriented "up" or "down", represented by the values +1 and -1. It has been used to study phase transitions in magnetic systems, where the system transitions from a low-temperature ordered state to a high-temperature disordered state. The Ising model has also been used to describe a wide range of phenomena, including order in biological systems, protein aggregation, neural network formation and also game theory. The energy of the system depends on the interactions between the spins and the energies of each particles, such that the Hamiltonian is:
$$
\widehat{H} = \lambda \sum_{j=1}^N \sigma^{(j)}_z + \sum_{j=1}^{N-1} \sigma^{(j)}_x \sigma^{(j+1)}_x,
$$
where the $\sigma_x$ and the $\sigma_z$ are the Pauli Matrix $$
\sigma_x = \begin{pmatrix}
0 & 1 \\
1 & 0
\end{pmatrix} \qquad \qquad \sigma_z = \begin{pmatrix}
1 & 0\\
0 & -1
\end{pmatrix}$$
The eigenvalues in the Ising model are distributed differently depending on the value of the parameter $\lambda$, which represents the strength of the interaction between the spins. In particular, for $\lambda$ values close to 0, the eigenvalues become degenerate, i.e. they cluster into a limited number of values, whereas for $\lambda$ values larger than 0, the eigenvalues begin to distinguish themselves.

The energy of the system is given by the sum of the interaction energies between the spins. In particular, each pair of neighboring spins contributes to the total energy of the system based on their interaction, which can be ferromagnetic or antiferromagnetic. In the first case, the spins tend to align parallelly, while in the second case they tend to align antiparallelly.

The Ising model can be numerically solved using the Hamiltonian. Solving the Hamiltonian allows the calculation of system properties, such as magnetization and free energy. Given a such Hamiltonian operator, the $j-$th particle contributes to that operator as $$
\widehat{H}_j  = \lambda \left( \mathbb{1}_{2\times 2}^{(1)} \otimes \mathbb{1}_{2 \times 2}^{(2)} \otimes \ldots \otimes \sigma_z^{(j)} \otimes \ldots \otimes \mathbb{1}_{2\times 2}^{(N)}\right) + \left( \mathbb{1}_{2 \times 2}^{(1)} \otimes \mathbb{1}_{2\times 2}^{(2)} \otimes \ldots \otimes \sigma_x^{(j)} \otimes \sigma_{x}^{(j+1)} \otimes \ldots \otimes \mathbb{1}_{2\times 2}^{(N)}\right). $$

This code numerically solves the Ising model by calculating the Hamiltonian of the system. The Hamiltonian is represented as a 2^N x 2^N matrix, where N is the number of particles, and is calculated separately for the interactions between the spins along the z-axis and along the x-axis. Subsequently, the eigenvalues and eigenvectors of the Hamiltonian matrix are calculated using the numpy.linalg.eigh function. The code also includes a progress bar to monitor the progress of the calculation and an estimate of the expected completion time.

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
