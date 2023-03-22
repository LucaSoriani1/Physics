import numpy as np
import math

# ======================================================================!
# This program aims to compute the algorithm on the renormalization
# group, namely ''real space renormalization group''or RSRG. The program 
# use this algorithm on a tranverse field Ising Model as exampel of quantum
# many-body problem. The program consider a fixed number of particles,
# then compute the Ising model. 
# ======================================================================!

# Hilbert space dimension
dimension = 2

# Size of the system = # of particles
N_particles = 2

# Bins for the strenght (lambda) parameter
binlambda = 300

# Precision for the evaluation of the results
precision=1e-4

# Max value of parameter lambda for the variation of the strenght parameter
p_lambda=3

# Discretization of the strenght parameter
dlambda = p_lambda/binlambda

def pauli():
    """
    It returns two matrices, one of which is the Pauli-X matrix and the other is the Pauli-Z matrix
    :return: Two matrices, one of which is the Pauli-X matrix and the other is the Pauli-Z matrix.
    """
    return [[0,1],[1,0]], [[1,0], [0,-1]]

def z_hamiltonian(sigma_z, N):    
    """
    The function takes in the number of particles and the sigma_z matrix and returns the 
    z-part of the Hamiltonian matrix

    :param sigma_z: the Pauli Z matrix
    :param N: number of particles
    :return: The z_hamiltonian is being returned.
    """

    # Creating the z-part of the Hamiltonian matrix of size 2^N x 2^N.
    z_ham = np.zeros((2**N, 2**N))
    
    # Starting the cycle about the computation of the Z-part.
    # I treat three cases: first when the index is i==0, i.e. I am
    # considering the first particle, the second when the index is
    # i==1, i.e. I am considering the second particle, and then
    # when the index is i>1. At the end of each cycle, The program
    # sums the result to the previous one. 
    for i in range(N):
        temp = np.eye(1)

        for j in range(i):
                temp = np.kron(temp, np.eye(2))

        temp = np.kron(temp, sigma_z)

        for j in range(N-i-1):
            temp = np.kron(temp, np.eye(2))        

        z_ham += temp

    return z_ham
        
def x_hamiltonian(sigma_x, N):
    """
    It creates the x-part for the Isign model that is the sum of all the possible combinations of the Pauli X matrix with the Identity.
    
    :param sigma_x: the Pauli X matrix
    :param N: number of particles in the system
    :return: the x-Hamiltonian for a given number of particles.
    """
    
    # Initialize the x-part of the Hamiltonian
    x_ham = np.zeros((2**N, 2**N))

    # Creating a Hamiltonian for a system of N.
    for i in range(N-1):

        # Creating a 1x1 matrix with 1s on the diagonal and 0s elsewhere.
        temp = np.eye(1)

        # Creating a tensor product of the identity matrix with itself.
        for j in range(i):
            temp = np.kron(temp, np.eye(2))

        # Double tensor product of the 'temp' with the sigma_x
        temp = np.kron(np.kron(temp, sigma_x), sigma_x)

        # The remianing tensor product with the identity
        for j in range(N-i-2):
            temp = np.kron(temp, np.eye(2))

        # Summing all the products
        x_ham += temp
        
    return x_ham

def projection(projector, hamiltonian):
    """
    It takes a projector and a Hamiltonian and returns the projected Hamiltonian
    
    :param projector: the projector matrix
    :param hamiltonian: the Hamiltonian of the system
    :return: The projection of the Hamiltonian onto the projector.
    """        
    return np.matmul(projector.T.conj(), np.matmul(hamiltonian, projector))


def reduced(rho, indk, dim):
    """
    This function computes the reduced density matrix of a given density matrix.
    The reduced density matrix is obtained by tracing out the degrees of freedom
    of the system that are not considered.
    """
    
    #Initialization of the reduced matrix as null matrix
    rows = int(math.sqrt(len(rho)))
    columns = int(math.sqrt(len(rho[0])))
    redu = np.zeros((rows, columns))
    
    #Compute the reduced density matrix
    for ii in range(rows):
        for jj in range(columns):
            for kk in range(dim):
                #row
                row = 1 + (ii % (dim**(indk-1))) + (kk-1)*(dim*(indk-1)) + dim*(ii-(ii % (dim**(indk-1)))) + (kk-1)*(2-indk)

                #col
                col = 1 + (jj % (dim**(indk-1))) + (kk-1)*(dim*(indk-1)) + dim*(jj-(jj % (dim**(indk-1)))) + (kk-1)*(2-indk)

                #reduced
                redu[ii,jj] = redu[ii,jj] + rho[row-1, col-1]
    return redu

def expected_time(times, it, cicles):
    """
    It calculates the expected time to complete a given number of iterations, based on the average time
    of the previous iterations
    
    :param times: the array of times for each iteration
    :param it: number of iterations
    :param cicles: number of cicles to run
    :return: The expected time to complete the remaining cicles.
    """
    sum = np.sum(times)

    return (((sum/it)*cicles)/60 - sum/60)

       

def progress_bar(progress, total, times, label=None):
    """
    It prints a progress bar with a remaining time estimate
    
    :param progress: the current progress
    :param total: the total number of iterations
    :param times: The time it took to complete the previous iteration
    :param label: The label to display before the progress bar
    """

    # Converting the time in minutes to seconds.
    minutes = int(times)
    seconds = int((times - minutes)*60)

    # Defining the percentage of progress
    percent = 100*(progress/float(total))
    
    # The above code is creating a variable called bar and assigning it a value of a string of '█'
    # multiplied by the integer of the percent variable. Then it is adding a string of '-' multiplied
    # by the integer of 100 minus the percent variable.
    bar = '█' * int(percent) + '-' * (100-int(percent))
    
    # Printing the progress bar and the remaining time.
    if progress > 1:
        if minutes > 0 :
            print(f"\033[1A\033[2K\rRemaining time: {minutes} minutes and {seconds} seconds\n|{bar}| {percent:.2f}% ({progress}/{total})", end="\r")
        else:
            print(f"\033[1A\033[2K\rRemaining time: {seconds} seconds\n|{bar}| {percent:.2f}% ({progress}/{total})", end="\r")

    else:
        if minutes >0 :
            print(f"Remaining time: {minutes} minutes and {seconds} seconds\n|{bar}| {percent:.2f}% ({progress}/{total})", end="\r")
        else:
            print(f"Remaining time: {seconds} seconds\n|{bar}| {percent:.2f}% ({progress}/{total})", end="\r")
            
    if progress == total:
        print(f"\033[1A\033[2K\rComplete!\n|{bar}| {percent:.2f}% ({progress}/{total})", end="\r")
        print("")