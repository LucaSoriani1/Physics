# Bose-Hubbard Model

The Bose-Hubbard model is a theoretical model used in condensed matter physics to describe the behavior of ultracold atoms trapped in optical lattices. The model was proposed by Satyendra Nath Bose and Jagdish Chandra Bose in the early 1920s, and was later generalized and applied to optical lattices by P.W. Anderson in the 1980s.

The Bose-Hubbard model describes a lattice of sites, each of which can be occupied by zero, one, or multiple identical bosonic atoms. The model considers two types of energy interactions: the kinetic energy of the atoms and the interactions between the atoms at neighboring sites. The kinetic energy is represented by the hopping parameter, which determines the probability of an atom moving from one site to a neighboring site. The interaction energy is represented by the on-site interaction parameter, which determines the energy cost or gain of having multiple atoms occupying the same site.

The Bose-Hubbard model is described by the Hamiltonian: $$
\widehat{H}_{B-H} = J\left(\widehat{b}_1^{\, \dagger\,}\widehat{b}_2 + \widehat{b}_1 \widehat{b}_2^{\,\dagger}\right) + \frac{U}{2}\left(\widehat{n}_1(\widehat{n}_1-\mathbb{I}) + \widehat{n}_2(\widehat{n}_2-\mathbb{I})\right) - \mu(\widehat{n}_1 + \widehat{n}_2),$$
where $\widehat{b}_i$ and $\widehat{b}^{\dagger}_i$ are the bosonic creation and annihilation operators at site $i$, $\widehat{n}_i$ are the particles-number operator for each site, $J$ is the hopping parameter, $U$ is the on-site interaction parameter, and $\mu$ is the chemical potential. The first term describes the hopping of atoms between neighboring sites, the second term describes the on-site interactions between atoms, and the third term describes the energy cost of having atoms present at each site.

The Bose-Hubbard model has been used to study a wide range of phenomena, including superfluidity, Mott insulator phases, and the Bose glass phase. It has also been used as a theoretical framework for designing and analyzing experiments involving ultracold atoms trapped in optical lattices.

Overall, the Bose-Hubbard model is a powerful tool for understanding the behavior of ultracold atoms in optical lattices, and has important applications in both fundamental physics and quantum technologies.

The provided code defines a function BoseHubbardHamiltonian that takes three parameters $J$, $U$, and $\mu$ and returns the Bose-Hubbard Hamiltonian for two bosons in a 2-site chain. The Hamiltonian is defined in terms of the hopping parameter, on-site interaction strength, and chemical potential.

The function calculates the Hamiltonian using the ladder operators, particle-number operator, and the creation and annihilation operators for each site. The function eigen_calculator takes the ground state as an input and computes the standard deviation of the number of particles for each site using the particle-number operator, such as $$
\left(\Delta n_1 \right)^2 = \bra{\psi_{gs}} (\widehat{n}_1)^2 \ket{\psi_{gs}} - (\bra{\psi_{gs}} \widehat{n}_1 \ket{\psi_{gs}})^2 = \langle \widehat{n}_1^2 \rangle_{gf} - (\langle \widehat{n}_1 \rangle_{gs})^2.$$

The function print_details prints the ground state in a user-friendly way. It takes three parameters: gs_index, gs, and delta. The function uses a dictionary to map the indices of the eigenstates to their corresponding quantum states, and then constructs a string to print the ground state and the amplitude of probability for each state. The function also prints the total probability of finding the bosons in the ground state.

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