# Real Space Renormalization Group (RSRG)
The RSRG algorithm consists of doubling the number of particles at each step. Once I get the doubled system, the computation of the ground state of the doubled system. After that, the matrix of the doubled system is projected to the system of the $N$ initial size, but it will contain the information of the doubled system. After that, the algorithm starts again: it doubles the system and computes the ground state, and then the system is projected again. The algorithm stops when the ground state is stable, thus when the last ground state is similar to the previous one. The steps for this algorithm are
1. The algorithm stars when the Hamiltonian operator is defined by the Ising Model, i.e. $H_{N} = \lambda\sum_{i=1}^N \sigma^{(i)}_z + \sum_{i=1}^{N-1} \sigma_x^{(i)} \sigma_x^{(i+1)}$. Now the algorithm doubles the system, introducing another Hamiltonian operator defined by $$
H_{2N} =
 H_{N} \otimes \mathbb{1}_{N} + \mathbb{1}_N \otimes H_N + 
\underbrace{ \mathbb{1}_{N-1} \otimes \sigma_x^{(N)} \otimes \sigma_x^{(N+1)} \otimes \mathbb{1}_{N-1}}_{H_{2N}^{int}}
$$
where the doubled system is formed by two specular systems near each other and the $H_{2N}^{int}$ is the interaction between the N-th particle of the first system and the first particle of the second system. In this way, the matrix that represents the Hamiltonian will have the dimensions of $2^{2N} \times 2^{2N}$, while the matrix of the initial Hamiltonian $H_{N}$ has dimension $2^N \times 2^N$. Thus a system of $N=2$ particles will produce an initial matrix $H_2$ of dimensions $4\times 4$ and the double system of $N=4$ particles will produce a matrix $H_4$ of dimensions $16 \times 16$.
2. Thus, once I get the doubled system, the program can compute the eigenvalues. The ground state energy is the quantity that I will evaluate. While the size of the system increases, the ground state decreases thus, in order to obtain the plateau for the ground state I must divide the eigenvalue for the number of particles for that cycle. More precisely $E_{g.s.}^k = E_0^k/N$ for each cycle $k$. 
3. The algorithm ends when the last ground state energy remains constat, thus the condition $|E^k_{g.s.} - E^{k-1}_{g.s.}| < \delta$ , where $\delta$ is an arbitrary tolerance.
4. If the condition in $3.$ is not respected, the algorithm proceeds computing the projector $P$: this is a matrix of dimensions $2^{2N} \times 2^N$ formed by the first $2^N$ eigenvectors. Once obtained $P$, the algorithm computes the projected Hamiltonian such that $$
\widetilde{H}_{N} = P^{\dagger} H_{2N} P = P^\dagger \left( H_{N} \otimes \mathbb{1}_{N} + \mathbb{1}_N \otimes H_N + 
\mathbb{1}_{N-1} \otimes \sigma_x^{(N)} \otimes \sigma_x^{(N+1)} \otimes \mathbb{1}_{N-1}\right)P.
$$
In this case, the system is projected to a sub-system that contains the information of the doubled one, thus it is written on another basis. The algorithm restarts from the point $1.$ of the algorithm, where now there is $\widetilde{H}_N$ instead of $H_N$. The doubled Hamiltonian will be $$
\widetilde{H}_{2N} = \widetilde{H}_N \otimes \mathbb{1}_{N} + \mathbb{1}_N \otimes \widetilde{H}_N + \mathbb{1}_{N-1} \otimes \widetilde{\sigma}_{x}^{(N)} \otimes \widetilde{\sigma}_{x}^{(N+1)} \otimes \mathbb{1}_{N-1}, $$
where I have to project the left and the right partof the $H_{2N}^{int}$ in order to obtain the term $\mathbb{1}_{N-1} \otimes \widetilde{\sigma}_{x}^{(N)}$ for the left part and the term $\widetilde{\sigma}_{x}^{(N+1)} \otimes \mathbb{1}_{N-1}$ for the right part written in the correct basis. This projection is made by $$
\mathbb{1}_{N-1} \otimes \widetilde{\sigma}_{x}^{(N)} = P^\dagger \left(\mathbb{1}_{N} \otimes \mathbb{1}_{N-1} \otimes \sigma_{x}^{(N)}\right)P \qquad \text{for the left part} \\
\widetilde{\sigma}_{x}^{(N+1)} \otimes \mathbb{1}_{N-1} = P^\dagger \left(\sigma_{x}^{(N+1)}\otimes \mathbb{1}_{N-1} \otimes \mathbb{1}_{N} \right)P \qquad \text{for the right part}$$
This algorithm will be repeated until the condition in $3.$ is true, then the algorithm stops because the system achieves the thermodynamic limit. 

# Infinite Density Matrix Renormalization Group (Infinite DMRG)

The infinite DMRG is another algorithm that aims to achieve the thermodynamic limit by adding just 2 particles between the two specular systems. The initial situation is analogous to the RSRG algorithm, i.e. the system is made by the Ising Model for $N$ particles. At this point, the algorithm is made to doubled the initial system and adding two particles between the two initial system. Once the program builds the total Hamiltonian for the system, it computes the diagonalization of the total Hamiltonian in order to check if the ground state is in the thermodynamic limit. If it is not, is needed to build the density matrix, starting from the eigenvector of the ground state. Now the algorithm requires to compute the reduced density matrix for the left and the right sub-system. These two reduced density matrices are used to build the projector, starting from their diagonalization, to project the left and right sub-system into another subsystem that contains the information about the previous sub-system with $N+1$ particles. The steps of this algorithm are

1. Starting from the total Hamiltonian of the Ising model $H_N$ with $N$ particles, the first step is to build the total Hamiltonian for $2N+2$ particles, namely $H_{2N+2}$. This system could be seen as two Ising Models, among which there are two particles. If $H_{N}$ is the Hamiltonian of the Ising Model,  the non-interaction part is made by $$
H_{2N+2}^{non-int} = H^1_{2N+2} + H^2_{2N+2} + H_{2N+2}^3 + H^4_{2N+2}  =  \\ 
= H_N \otimes \mathbb{1}_{N+2} + \mathbb{1}_N \otimes \lambda\cdot\sigma_x^{(N+1)} \otimes \mathbb{1}_{N+1} + \mathbb{1}_{N+1} \otimes \lambda\cdot\sigma_z^{(N+2)} \otimes \mathbb{1}_{N} + \mathbb{1}_{N+2} \otimes H_N
$$
In the first and the last term, I recognize the initial Hamiltonian $H_N$ for the two separated systems. The second and the third term are contributions of the single particles at position $N+1$ and $N+2$, between the two Ising Model. These two terms are formed by the strength term of the Ising Model because they include only one particle per term. 
Now, the interaction part is made by $$
H_{2N+2}^{int} = H_{2N+2}^{1-2} + H_{2N+2}^{2-3} + H^{3-4}_{2N+2} =  \\ = \mathbb{1}_{N-1} \otimes \sigma_x^{(N)} \otimes \sigma_x^{(N+1)} \otimes \mathbb{1}_{N+1} +\mathbb{1}_{N} \otimes \sigma_x^{(N+1)} \otimes \sigma_x^{(N+2)} \otimes \mathbb{1}_{N} +\mathbb{1}_{N+1} \otimes \sigma_x^{(N+2)} \otimes \sigma_x^{(N+3)} \otimes \mathbb{1}_{N-1}$$
In the end, the interaction term is the interaction between the $N$-th particle of the first system and the first free particle, the interaction between the two free particles, and the interaction between the second free particle and the first of the second system. 
Now, the program can build the $H_{2N+2}$ of the total system as $H_{2N+2} = H_{2N+2}^{non-int} + H_{2N+2}^{int}$.

2. Now there is the diagonalization of the $H_{2N+2}$ Hamiltonian. The first eigenvalue divided by the number of particles will be the ground state of the system, namely $E^k_{g.s.} = E^k_0/N$, where $k$ is the cycle. 

3. As seen for the point $3.$ in the RSRG algorithm, the program will prove is the ground state energy is constant with the condition $|E^k_{g.s.} - E^{k-1}_{g.s.}|<\delta$, where $\delta$ is an arbitrary tolerance. 

4. if the previous condition is not respected, the program takes the eigenvector associated to the ground state energy and it builds the density matrix $\rho = \ket{\psi_{g.s.}} \bra{\psi_{g.s.}}$.

5. From the density matrix $\rho$, the program computes the reduced density matrix for the left-system and the right system. For example, for the left system, the reduced density matrix is given by $$
\rho_L = Tr_R \,  \rho $$
and the reduced density matrix is diagonalized.

6. From the eigenvectors of the diagonalized reduced density matrix, the program builds the projector $P$. This projector is a matrix of size $2^{N+1} \times 2^{N}$ and it is formed by the eigenvectors of the reduced density matrix stored in the descending order.

7. The projector $P$ projects the temrs $H^L_{N+1} = H^1_{N+1} + H^2_{N+1} + H^{1-2}_{N+1}$ and the term $L_{N+1} = \mathbb{1}_{N} \otimes \sigma_x^{(N+1)}$ into the sub-system of $N$ particles, but with the information of the system formed by $N+1$ particels. In analogous, the right part formed by the $H^R_{N+1} = H^3_{N+1} + H^4_{N+1} + H^{3-4}_{N+1}$  and the $R_{N+1} = \sigma_x^{(N+2)} \otimes \mathbb{1}_{N}$ is procjeted. And then the algorithm restarts from the point $1)$ with the new Hamiltonians given by $$
\widetilde{H}^L_{N} = P^\dagger \left(H^L_{N+1}\right) P 
\qquad \widetilde{H}^R_{N} = P^\dagger \left(H^R_{N+1}\right) P  \\    \\
\mathbb{1}_{N-1} \otimes \widetilde{\sigma}_x^{(N)} = P^\dagger L_{N+1} P \qquad \widetilde{\sigma}_x^{(N+1)} \otimes \mathbb{1}_{N-1} = P^\dagger R_{N+1} P $$
The algorithm will be repeated until the condition at the point $3.$ is verified. At this point, the thermodynamic limit is achieved and the program stops the computation.

### Code implementation

Renormalization group theory (RG) is a physical theory that describes how the properties of a physical system change when the system is examined at different length or energy scales. In particular, RG is used to study physical systems that have a multi-scale structure, such as critical systems, turbulent fluids and disordered materials. 

Renormalization group theory was initially developed for statistical physics, but was later extended to other areas of physics, such as particle physics, string theory and condensed matter physics. 

This theory involves the construction of a group of transformations that change the length or energy scale of the physical system. The properties of the system at one scale are then described by the properties of the system at a different scale. This process is repeated iteratively until the desired scale is reached. 

An important aspect is the notion of a fixed point. A fixed point is a set of conditions in which the properties of the system remain unchanged under the transformations of the renormalization group. The system evolves autonomously around a fixed point, and its properties encrypt information about the physics of the system at all scales. 

Real space renormalization group (RSRG) and infinite density matrix renormalization group (Infinite DMRG) are two computational techniques that use renormalization group theory to study many-body quantum systems. In this example we are going to apply these two techniques to obtain the "fixed point" ground state, looking at efficiency and accuracy. As a many-body system, we will take the Isign model. 

## RSRG

The RSRG algorithm consists of fixing a number of particles and doubling them at each iteration. Once the doubled system is obtained, the fundamental state of the doubled system is calculated. Next, the matrix of the doubled system is projected onto the system of initial dimension $N$, but it will contain the information of the doubled system. After that, the algorithm starts again: it doubles the system and calculates the fundamental state, then the system is projected again. The algorithm stops when the fundamental state is stable, so when the last fundamental state is similar to the previous one. The steps of this algorithm are:

1. The algorithm stars when the Hamiltonian operator is defined by the Ising Model, i.e. $H_{N} = \lambda\sum_{i=1}^N \sigma^{(i)}_z + \sum_{i=1}^{N-1} \sigma_x^{(i)} \sigma_x^{(i+1)}$. Now the algorithm doubles the system, introducing another Hamiltonian operator defined by 
$$H_{2N} =
 H_{N} \otimes \mathbb{1}_{N} 
 + \mathbb{1}_N \otimes H_N + 
\underbrace{ \mathbb{1}_{N-1} \otimes \sigma_x^{(N)} \otimes \sigma_x^{(N+1)} \otimes \mathbb{1}_{N-1}}_{H_{2N}^{int}}$$

where the doubled system is formed by two specular systems near each other and the $H_{2N}^{int}$ is the interaction between the N-th particle of the first system and the first particle of the second system. In this way, the matrix that represents the Hamiltonian will have the dimensions of $2^{2N} \times 2^{2N}$, while the matrix of the initial Hamiltonian $H_{N}$ has dimension $2^N \times 2^N$. Thus a system of $N=2$ particles will produce an initial matrix $H_2$ of dimensions $4\times 4$ and the double system of $N=4$ particles will produce a matrix $H_4$ of dimensions $16 \times 16$.
2. Thus, once I get the doubled system, the program can compute the eigenvalues. The ground state energy is the quantity that I will evaluate. While the size of the system increases, the ground state decreases thus, in order to obtain the plateau for the ground state I must divide the eigenvalue for the number of particles for that cycle. More precisely $E_{g.s.}^k = E_0^k/N$ for each cycle $k$.
3. The algorithm ends when the ground state energy remains constat, thus the condition $|E^k_{g.s.} - E^{k-1}_{g.s.}| < \delta$ , where $\delta$ is an arbitrary tolerance.
4. If the condition in $3.$ is not respected, the algorithm proceeds computing the projector $P$: this is a matrix of dimensions $2^{2N} \times 2^N$ formed by the first $2^N$ eigenvectors. Once obtained $P$, the algorithm computes the projected Hamiltonian such that

$$\widetilde{H}_{N} = P^{\dagger} H_{2N} P = P^\dagger \left( H_{N} \otimes \mathbb{1}_{N} 
 + \mathbb{1}_N \otimes H_N + 
\mathbb{1}_{N-1} \otimes \sigma_x^{(N)} \otimes \sigma_x^{(N+1)} \otimes \mathbb{1}_{N-1}\right)P.$$

In this case, the system is projected to a sub-system that contains the information of the doubled one, thus it is written on another basis. The algorithm restarts from the point $1.$ of the algorithm, where now there is $\widetilde{H}_N$ instead of $H_N$. The doubled Hamiltonian will be

$$
\widetilde{H}_{2N} = \widetilde{H}_N \otimes \mathbb{1}_{N} + \mathbb{1}_N \otimes \widetilde{H}_N + \mathbb{1}_{N-1} \otimes \widetilde{\sigma}_{x}^{(N)} \otimes \widetilde{\sigma}_{x}^{(N+1)} \otimes \mathbb{1}_{N-1},$$

where it needs to project the left and the right part of the $H_{2N}^{int}$ in order to obtain the term $\mathbb{1}_{N-1} \otimes \widetilde{\sigma}_{x}^{(N)}$ for the left part and the term $\widetilde{\sigma}_{x}^{(N+1)} \otimes \mathbb{1}_{N-1}$ for the right part written in the correct basis. This projection is made by 

$$
\mathbb{1}_{N-1} \otimes \widetilde{\sigma}_{x}^{(N)} = P^\dagger \left(\mathbb{1}_{N} \otimes \mathbb{1}_{N-1} \otimes \sigma_{x}^{(N)}\right)P $$

 for the left part, while
$$
\widetilde{\sigma}_{x}^{(N+1)} \otimes \mathbb{1}_{N-1} = P^\dagger \left(\sigma_{x}^{(N+1)}\otimes \mathbb{1}_{N-1} \otimes \mathbb{1}_{N} \right)P$$

This algorithm will be repeated until the condition in $3.$ is true, then the algorithm stops because the system achieves the thermodynamic limit.

## Infinite DMRG

The infinite DMRG is another algorithm that aims to achieve the thermodynamic limit by adding just $2$ particles between the two specular systems. The initial situation is analogous to the RSRG algorithm, i.e. the system is made by the Ising Model for $N$ particles. At this point, the algorithm is made to doubled the initial system and adding two particles between the two initial system. Once the program builds the total Hamiltonian for the system, it computes the diagonalization of the total Hamiltonian in order to check if the ground state is in the thermodynamic limit. If it is not, is needed to build the density matrix, starting from the eigenvector of the ground state. Now the algorithm requires to compute the reduced density matrix for the left and the right sub-system. These two reduced density matrices are used to build the projector, starting from their diagonalization, to project the left and right sub-system into another subsystem that contains the information about the previous sub-system with $N+1$ particles. The steps of this algorithm are


1. Starting from the total Hamiltonian of the Ising model $H_N$ with $N$ particles, the first step is to build the total Hamiltonian for $2N+2$ particles, namely $H_{2N+2}$. This system could be seen as two Ising Models, among which there are two particles. If $H_{N}$ is the Hamiltonian of the Ising Model,  the non-interaction part is made by 
$$ H_{2N+2}^{non-int}  H^1_{2N+2} + H^2_{2N+2} + H_{2N+2}^3 + H^4_{2N+2}   
= H_N \otimes \mathbb{1}_{N+2} + \mathbb{1}_N \otimes \lambda\cdot\sigma_x^{(N+1)} \otimes \mathbb{1}_{N+1} + \mathbb{1}_{N+1} \otimes \lambda\cdot\sigma_z^{(N+2)} \otimes \mathbb{1}_{N} + \mathbb{1}_{N+2} \otimes H_N
$$

In the first and the last term, it recognizes the initial Hamiltonian $H_N$ for the two separated systems. The second and the third term are contributions of the single particles at position $N+1$ and $N+2$, between the two Ising Model. These two terms are formed by the strength term of the Ising Model because they include only one particle per term. 
Now, the interaction part is made by 
$$ H_{2N+2}^{int} = H_{2N+2}^{1-2} + H_{2N+2}^{2-3} + H^{3-4}_{2N+2}  = \mathbb{1}_{N-1} \otimes \sigma_x^{(N)} \otimes \sigma_x^{(N+1)} \otimes \mathbb{1}_{N+1} +\mathbb{1}_{N} \otimes \sigma_x^{(N+1)} \otimes \sigma_x^{(N+2)} \otimes \mathbb{1}_{N} +\mathbb{1}_{N+1} \otimes \sigma_x^{(N+2)} \otimes \sigma_x^{(N+3)} \otimes \mathbb{1}_{N-1}

In the end, the interaction term is the interaction between the $N$-th particle of the first system and the first free particle, the interaction between the two free particles, and the interaction between the second free particle and the first of the second system. 

Now, the program can build the $H_{2N+2}$ of the total system as $H_{2N+2} = H_{2N+2}^{non-int} + H_{2N+2}^{int}$.
2. Now there is the diagonalization of the $H_{2N+2}$ Hamiltonian. The first eigenvalue divided by the number of particles will be the ground state of the system, namely $E^k_{g.s.} = E^k_0/N$, where $k$ is the cycle. 
3. As seen for the point $3.$ in the RSRG algorithm, the program will prove is the ground state energy is constant with the condition $|E^k_{g.s.} - E^{k-1}_{g.s.}|<\delta$, where $\delta$ is an arbitrary tolerance. 
4. If the previous condition is not respected, the program takes the eigenvector associated to the ground state energy and it builds the density matrix $\rho = \ket{\psi_{g.s.}} \bra{\psi_{g.s.}}$.
5. From the density matrix $\rho$, the program computes the reduced density matrix for the left-system and the right system. For example, for the left system, the reduced density matrix is given by $$
\rho_L = Tr_R \,  \rho
$$
and the reduced density matrix is diagonalized. 
6. From the eigenvectors of the diagonalized reduced density matrix, the program builds the projector $P$. This projector is a matrix of size $2^{N+1} \times 2^{N}$ and it is formed by the eigenvectors of the reduced density matrix stored in the descending order.
7. The projector $P$ projects the terms $H^L_{N+1} = H^1_{N+1} + H^2_{N+1} + H^{1-2}_{N+1}$ and the term $L_{N+1} = \mathbb{1}_{N} \otimes \sigma_x^{(N+1)}$ into the sub-system of $N$ particles, but with the information of the system formed by $N+1$ particels. In analogous, the right part formed by the $H^R_{N+1} = H^3_{N+1} + H^4_{N+1} + H^{3-4}_{N+1}$  and the $R_{N+1} = \sigma_x^{(N+2)} \otimes \mathbb{1}_{N}$ is procjeted. And then the algorithm restarts from the point $1)$ with the new Hamiltonians given by 

$$
\widetilde{H}^L_{N} = P^\dagger \left(H^L_{N+1}\right) P 
\qquad \quad  \widetilde{H}^R_{N} = P^\dagger \left(H^R_{N+1}\right) P  $$

$$
\mathbb{1}_{N-1} \otimes \widetilde{\sigma}_x^{(N)} = P^\dagger L_{N+1} P \qquad \quad \widetilde{\sigma}_x^{(N+1)} \otimes \mathbb{1}_{N-1} = P^\dagger R_{N+1} P
$$

The algorithm will be repeated until the condition at the point $3.$ is verified. At this point, the thermodynamic limit is achieved and the program stops the computation.


The implementation of the script is based on the numpy library, and also imports several utility functions from the 'renormalization_utils' module. This module provides all the tools to develop the algorithms, such as creating the Hamiltonian of the Isign model and the function responsible for projection. The renormalization algorithm is used to compute the ground state energy of the Ising model as a function of the coupling strength parameter $\lambda$.

The function first initializes several parameters:

* dimension: the Hilbert space dimension of each particle, default is $2$.
* N_particles: the size of the system (number of particles), default is $2$.
* binlambda: the number of bins for the strength parameter$\lambda$, default is $200$.
* precision: the precision for the evaluation of the results, default is $10^{-5}$.
* p_lambda: the maximum value of the strength parameter$\lambda$, default is $3$.

The function then initializes the Pauli matrices, which are used to define the interaction and non-interaction parts of the Ising model. The function then initializes several arrays to store the partial time of the computation, the steps and levels at $\lambda=0$ and $\lambda=max$, the convergence steps, $\lambda$, and ground state energy.

The function then enters a loop over the bins for the strength parameter. Within each iteration of the loop, the function computes the full Hamiltonian of the Ising model and initializes the left and right interaction terms. The function then enters a while loop to perform the renormalization algorithm. The function then computes the eigenvalues and eigenvectors of the doubled Hamiltonian, and stores the ground state energy. The function then checks for convergence by comparing the ground state energy of the current iteration with the ground state energy of the previous iteration. If the ground state energies are similar, the function stops iterating and stores the final results. If the ground state energies are not similar, the function updates the full Hamiltonian, left interaction term, and right interaction term by projecting the doubled Hamiltonian onto the previous iteration's basis. This projection is done using the utility functions from the renormalization_utils module. The function then stores the new $\lambda$ and ground state energy, and increments the step counter.

After the while loop completes, the function stores the final results for the current bin of $\lambda$, and moves to the next bin. Once all the bins have been iterated over, the function returns the results, which consist of the steps and levels at $\lambda=0$ and $\lambda=max$, the convergence steps, $\lambda$, and ground state energy. These results can then be used to analyze the behavior of the Ising model under renormalization.

### Usage
To run the program, open the terminal and type the following command:
```
python random_trajectories.py
```
Five random trajectories will be generated and plotted.

### Dependencies
The program requires the installation of the following Python libraries:

* numpy
* matplotlib
* scipy

The libraries can be installed using pip, for example:
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
