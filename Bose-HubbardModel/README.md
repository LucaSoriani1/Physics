# Bose-Hubbard Model

The Bose-Hubbard model is a theoretical model used in condensed matter physics to describe the behavior of ultracold atoms trapped in optical lattices. The model was proposed by originally by John Hubbard, while the term 'Bose' refers to the fact that this particular application is composed of bosonic particles.

The Bose-Hubbard model describes a lattice of potential holes, each of which can be occupied by zero, one or more identical bosonic atoms. The model considers two types of energy interactions: the kinetic energy of atoms and the interactions between atoms at the same site. Kinetic energy is represented by the hopping parameter, which determines the probability of an atom moving from one site to a nearby site. Interaction energy is represented by the in-site interaction parameter, which determines the energy cost or gain of having multiple atoms occupying the same site.

The Bose-Hubbard model is described by the Hamiltonian: $$
\widehat{H}_{B-H} = J\left(\widehat{b}_1^{\, \dagger\,} \widehat{b}_2 + \widehat{b}_1 \widehat{b}_2^{\dagger}\right) + \frac{U}{2}\left(\widehat{n}_1(\widehat{n}_1-\mathbb{I}) + \widehat{n}_2(\widehat{n}_2-\mathbb{I})\right) - \mu(\widehat{n}_1 + \widehat{n}_2), $$

where $\widehat{b}_i$ and $\widehat{b}^{\dagger}_i$ are the bosonic creation and annihilation operators at site $i$, $\widehat{n}_i$ are the particle-number operators for each site, $J$ is the hopping parameter, $U$ is the particle interaction parameter for each site, and $\mu$ is the chemical potential affecting the total number of particles in the model. The first term describes the hopping of atoms between neighboring sites, the second term describes the in situ interactions between atoms, and the third term describes the energy cost of the presence of atoms.

The Bose-Hubbard model has been used to study a wide range of phenomena, including superfluidity, Mott isolation phases and the Bose glass phase. It has also been used as a theoretical framework to design and analyze experiments involving ultracold atoms trapped in optical lattices.

The different regimes of the model depend on the parameters, such as, when the strength of the interactions at site $U$ is much greater than the jump parameter $J$, the system is in a Mott isolation phase, in which the atoms are bound to individual sites and cannot move easily. Conversely, when the jump parameter $J$ is much greater than the in-site interaction parameter $U$, atoms can move easily between sites and the system is in a superfluid state.

	This model is of great importance in the realization of quantum technologies, such as quantum simulation and quantum computation. Experiments with ultracold atoms trapped in optical lattices are becoming increasingly precise and sophisticated, and the Bose-Hubbard model continues to be an important tool for understanding and developing these technologies.

The script in question was used as support for the calculations in the dissertation. Specifically, I assumed that I had two particles in two potential holes, so that I had the following operators: $$
\widehat{b}=\left(\begin{array}{ccc}{0} & {\sqrt{1}} & {0} \\ {0} & {0} & {\sqrt{2}} \\   {0} & {0} & {0}\end{array}\right),
\, \, \,
\widehat{b}^{\, \dagger}=\left(\begin{array}{ccc}{0} & {0} & {0} \\ {\sqrt{1}} & {0} & {0} \\ {\sqrt{2}} & {0} & {0}\end{array}\right), \, \, \,
\widehat{\mathcal{N}}=\left(\begin{array}{ccc}{0} & {0} & {0} \\ {0} & {1} & {0}  \\ {0} & {0} & {2} \end{array}\right).$$

The constructed model is aimed at demonstrating the phase transition from Mott insulator to superfluid state and vice versa. Calculations were performed by varying the relationship between the hopping parameter $J$ and the interaction potential $U$. Once the interaction $U=1$ is fixed for simplicity and having defined the chemical potential as $\mu = U/2 = 1/2$, the script goes on to calculate the groundstate of the system for different values of $J$. Theoretically, the groundstate of such a system can be composed of 9 different states: $ket{00}$, $ket{01}$, $ket{02}$, $ket{10}$, $ket{11}$, $ket{12}$, $ket{20}$, $ket{2,1}$, $ket{2,2}$, where the first value represents the number of particles in the first potential hole and the second the number of particles in the second. Physically, however, only 3 of them are possible, namely those whose sum of particles is equal to $2$, that is, the states $\ket{20}$, $\ket{11}$$ and \ket{02}$. Thus the assumption of the values of $J$, $U$ and $\mu$ was chosen to meet this constraint.

The phase change can be mathematically verified by looking at the fluctuation of the number of particles for each value of $J/U$, so that: $$
		
\left(\Delta n_1 \right)^2 = \bra{\psi_{gs}} (\widehat{n}_1)^2 \ket{\psi_{gs}} - (\bra{\psi_{gs}} \widehat{n}_1 \ket{\psi_{gs}})^2 = \langle \widehat{n}_1^2 \rangle_{gs} - (\langle \widehat{n}_1 \rangle_{gs})^2$$

To do this, the code defines a BoseHubbardHamiltonian(J, U, mu) function that takes as input three parameters, the hopping parameter J, the on-site interaction strength U, and the chemical potential mu, and returns the Bose-Hubbard Hamiltonian for a chain of 2 sites with 2 bosons. The Hamiltonian is defined in terms of the creation and destruction operators and the number-particle operator for each site.

Then the function eigen_calculator(gs) is defined, which takes the ground state as input and calculates the standard deviation of the number of particles for each site using the particle counting operator.


The code uses the numpy library for mathematical operations and the tabulate library for printing the table.

Finally, the function print_details(gs_index, gs, delta) is defined that prints the fundamental state in a user-friendly way. The function uses a dictionary to map the indices of the eigenstates to their corresponding quantum states and constructs a string to print the fundamental state and probability amplitude for each state. The function also prints the total probability of finding bosons in the fundamental state.

### Usage
To run the program, open the terminal and type the following command:
```
python hubbard.py
```

The default parameters are:
* U=1
* mu=U/2=1/2

and they are taken for semplicity in order to compute phase transition between Mott insulator and superfluid.

### Dependencies
The program requires the installation of the following Python libraries:

* numpy

The libraries can be installed using pip, for example:
```
pip install -r requirements.txt
```

### Modifications and Contributions
This program can be extended and modified to suit different purposes, such as modeling natural processes or generating random data for testing purposes. If you would like to contribute to the program or propose modifications, please send a pull request or contact the author.