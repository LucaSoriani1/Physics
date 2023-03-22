import numpy as np
from tabulate import tabulate

"""
The Bose-Hubbard model is a theoretical model used in condensed matter physics
to describe the behavior of ultracold atoms trapped in optical lattices.

The Bose-Hubbard model describes a lattice of sites,
each of which can be occupied by zero, one, or multiple identical bosonic atoms.
The model considers two types of energy interactions:
the kinetic energy of the atoms and the interactions between the atoms at neighboring sites.
The kinetic energy is represented by the hopping parameter,
which determines the probability of an atom moving from one site to a neighboring site.
The interaction energy is represented by the on-site interaction parameter,
which determines the energy cost or gain of having multiple atoms occupying the same site.

The Bose-Hubbard model is described by the Hamiltonian:

H = -J ∑(b†i b_{i+1} + b†_{i+1} b_i) + U/2 ∑ni(ni-I) - μ ∑ni,

where bi and b†i are the bosonic creation and annihilation operators at site i,
ni are the particles-number operator for each site,
J is the hopping parameter, U is the on-site interaction parameter, and μ is the chemical potential.
The first term describes the hopping of atoms between neighboring sites,
the second term describes the on-site interactions between atoms,
and the third term describes the energy cost of having atoms present at each site.

The Bose-Hubbard model has been used to study a wide range of phenomena,
including superfluidity, Mott insulator phases, and the Bose glass phase.
It has also been used as a theoretical framework for designing and
analyzing experiments involving ultracold atoms trapped in optical lattices.

Overall, the Bose-Hubbard model is a powerful tool for understanding
the behavior of ultracold atoms in optical lattices
and has important applications in both fundamental physics and quantum technologies.

The provided code defines a function BoseHubbardHamiltonian that takes three parameters J, U, and mu
and returns the Bose-Hubbard Hamiltonian for two bosons in a 2-site chain.
The Hamiltonian is defined in terms of the hopping parameter, on-site interaction strength, and chemical potential.

The function calculates the Hamiltonian using the ladder operators, particle-number operator,
and the creation and annihilation operators for each site.

The function eigen_calculator takes the ground state as an input
and computes the standard deviation of the number of particles for each site using the particle-number operator.

The function print_details prints the ground state in a user-friendly way.
It takes three parameters: gs_index, gs, and delta.
The function uses a dictionary to map the indices of the eigenstates to their corresponding quantum states,
and then constructs a string to print the ground state and the amplitude of probability for each state.
The function also prints the total probability of finding the bosons in the ground state.
"""


# Define the Bose-Hubbard Hamiltonian for a 2-site chain
def BoseHubbardHamiltonian(J, U, mu):
    """
    It takes the hopping parameter J, the on-site interaction parameter U, and the chemical
    potential mu as inputs, and returns the Bose-Hubbard Hamiltonian for 2 bosons in 2-site chain
    
    :param J: the hopping parameter
    :param U: the on-site interaction strength
    :param mu: chemical potential
    :return: The Hamiltonian for a 2-site chain.
    """
    # Define ladder operators
    b_create = np.array([[0,0,0],[1,0,0],[0,2**0.5,0]])
    b_destroy = np.array([[0,1,0],[0,0,2**0.5],[0,0,0]])

    #Define the particles-number operator
    n = np.matmul(b_create, b_destroy)

    # Define ladder operators for each site
    b_1_create = np.kron(b_create, np.eye(3))
    b_1_destroy = np.kron(b_destroy, np.eye(3))

    b_2_create = np.kron(np.eye(3), b_create)
    b_2_destroy = np.kron(np.eye(3), b_destroy)
    
    #Define the particles-number operator for each site
    n_1 = np.kron(n, np.eye(3))
    n_2 = np.kron(np.eye(3), n)
    
    
    #Define the hopping term for Hubbard Hamiltonian
    H_hop = J * (np.matmul(b_1_create, b_2_destroy) + np.matmul(b_1_destroy, b_2_create))

    #Define the confinement potential term for Hubbard Hamiltonian
    H_pot = (U/2) * (n_1*(n_1 - np.eye(9)) + n_2*(n_2 - np.eye(9)))


    #Define the chemical potential term for Hubbard Hamiltonian
    H_chem = mu * (n_1 + n_2)
    
    #Define the total Hamiltonian for Hubbard Hamiltonian
    H_hubbard = H_hop + H_pot + H_chem
    
    return H_hubbard

# From the eigenvalue and eigenvector, I find the ground state
# Computation of the Delta for the number of particles for each site
def eigen_calculator(gs):
    """
    The function takes the ground state as an input and returns the standard deviation of the number of
    particles for each site
    
    :param gs: the ground state
    :return: the standard deviation of the number of particles for each site.
    """
    
    # particles-number opterator
    n = np.array([[0,0,0], [0,1,0], [0,0,2]])
    
    # Delta n = (<gs|n^2|gs> - (<gs|n|gs>)^2)**0.5
    n_square = np.kron(np.matmul(n,n), np.eye(3))
    n_exp = np.kron(n, np.eye(3))
    
    delta_1 = np.dot(gs.T, n_square).dot(gs)
    delta_2 = np.dot(gs.T, n_exp).dot(gs) * np.dot(gs.T, n_exp).dot(gs)
    
    delta = (delta_1 - delta_2)**0.5
    
    return delta

def print_details(gs_index, gs, delta):
    """
    The function prints the ground state in a user-friendly way.
    It uses a dictionary to map the indices of the eigenstates to their corresponding quantum states,
    and then constructs a string to print the ground state and the amplitude of probability for each state.
    The function also prints the total probability of finding the bosons in the ground state.
    
    :param gs_index: the index of the ground state
    :param gs: the ground state
    :param delta: the standard deviation of the number of particles for each site
    """
    states={
        0:'|00>',
        1:'|01>',
        2:'|02>',
        3:'|10>',
        4:'|11>',
        5:'|12>',
        6:'|20>',
        7:'|21>',
        8:'|22>'
    }
    
    occuped_state = {}
    for i in range(len(gs)):
        if round(gs[i], 3) != 0:
            occuped_state[states[i]] = gs[i]
    
    table = []
    for k,v in occuped_state.items():
        table.append([k, round(v, 3)])
    ret = "\nEIGENVALUE FOR GROUND STATE: " + str(gs_index) + ',\n'
    ret += "GROUND STATE: "
    i=0
    sum=0
    for k,v in occuped_state.items():
        sum += v**2
        if v>0:
            if i!=0:
                ret += '+'
            ret+= str(v)[:5] + k
        else:
            ret += str(v)[:6] + k
            
        i+=1
    ret += ',\n\n'
    
    ret+=tabulate(table, headers=['State', 'Amplitude of Probability'])
    
    ret  += "\n\nTotal probability: " + str(sum)
    
    ret += "\n \nDelta of number of particles: " + str(delta) + "\n"
    print(ret)
    

def main():
    
    # Ask the user for the hopping parameter J
    J = float(input("Enter the hopping parameter J/U: "))

    # Set U = 1 for simplicity
    U = 1

    # Compute chemical potential from µ = U/2
    mu = U/2
    
    # Compute the Hubbard Hamiltonian
    H = BoseHubbardHamiltonian(-J, U, -mu)
    
    # Find the eigenvalues and eigenvectors    
    eigvals, eigvecs = np.linalg.eigh(H)
    
    # From the eigenvalue and eigenvector, I find the ground state
    # find the min eigenvalue and the eigenvector corrisponding to (ground state)
    gs_index = np.argmin(eigvals)
    gs = eigvecs[:, gs_index]
    
    delta = eigen_calculator(gs)
    
    print_details(eigvals[gs_index], gs, delta)
    
    
if __name__ == "__main__":
    main()

