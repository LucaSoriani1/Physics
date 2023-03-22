import matplotlib.pyplot as plt
import time
import numpy as np


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
# The program are made to compute the Hamiltonia operator for ''N''
# particles with spin 1/2 that form the Ising Model. In this model
# appears the Pauli's matrices and the interaction strength.
# The program computes the 2^N X 2^N representation for the Hamiltonian
# for different ''N'' and for the parameters between 0 and 3.
# Thus, the program compute the eigenvalues and eigenvectors problem.
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!


# The class Matrix() creates a square matrix of size dim x dim and initializes all elements and the trace to zero
class Matrix():

    def __init__(self, dim):
        """
        This function initializes a square matrix of size dim x dim and sets the trace to 0.

        :param dim: the dimension of the matrix
        """
        # Setting the row and column of the matrix to the same value.
        self.row = dim
        self.col = dim

        # Creating a matrix of size dim x dim and assigning the value of 0 to all the elements.
        self.matrix = np.zeros((dim, dim))

def pauli():
    """
    It creates two matrices, one for Pauli-x $\sigma_x$ and one for Pauli-z $\sigma_z$, and then calls the function
    `compute_trace` on each of them

    :return: The two Pauli matrices.
    """

    # Creating the Pauli-x 2x2 matrix and assigning the values of 1 to the elements in the first row and second
    # column and the second row and first column.
    sigma_x = Matrix(2)
    sigma_x.matrix[0,1] = 1
    sigma_x.matrix[1,0] = 1
    
    # Creating the Pauli-z 2x2 matrix and assigning the values of 1 to the element in the first row and first
    # column and -1 to the element in the second row and second column.
    sigma_z = Matrix(2)
    sigma_z.matrix[0,0] = 1
    sigma_z.matrix[1,1] = -1
    
    return sigma_x, sigma_z

# =====================================================================!
# This function creates the Z-part of the Ising Model. It takes as
# input the Sigma-Z Pauli's matrix and the integer 'n_particles',
# i.e. the number of particles (thus the iteraction).
# As output, the subroutine gives back the Z-part matrix.
# =====================================================================!
def z_hamiltonian(sigma_z, n_particles):    
    """
    The function takes in the number of particles and the sigma_z matrix and returns the 
    z-part of the Hamiltonian matrix

    :param sigma_z: the Pauli Z matrix
    :param n_particles: number of particles
    :return: The z_hamiltonian is being returned.
    """

    # Creating the z-part of the Hamiltonian matrix of size 2^n_particles x 2^n_particles.
    z_ham = Matrix(2**n_particles)
    
    # Starting the cycle about the computation of the Z-part.
    # I treat three cases: first when the index is i==0, i.e. I am
    # considering the first particle, the second when the index is
    # i==1, i.e. I am considering the second particle, and then
    # when the index is i>1. At the end of each cycle, The program
    # sums the result to the previous one. 
    for i in range(n_particles):
        temp = np.eye(1)

        for j in range(i):
                temp = np.kron(temp, np.eye(2))

        temp = np.kron(temp, sigma_z.matrix)

        for j in range(n_particles-i-1):
            temp = np.kron(temp, np.eye(2))        

        z_ham.matrix += temp

    return z_ham
        
def x_hamiltonian(sigma_x, n_particles):
    """
    It creates the x-part for the Isign model that is the sum of all the possible combinations of the Pauli X matrix with the Identity.
    
    :param sigma_x: the Pauli X matrix
    :param n_particles: number of particles in the system
    :return: the x-Hamiltonian for a given number of particles.
    """
    
    # Initialize the x-part of the Hamiltonian
    x_ham = Matrix(2**n_particles)

    # Creating a Hamiltonian for a system of n_particles.
    for i in range(n_particles-1):

        # Creating a 1x1 matrix with 1s on the diagonal and 0s elsewhere.
        temp = np.eye(1)

        # Creating a tensor product of the identity matrix with itself.
        for j in range(i):
            temp = np.kron(temp, np.eye(2))

        # Double tensor product of the 'temp' with the sigma_x
        temp = np.kron(np.kron(temp, sigma_x.matrix), sigma_x.matrix)

        # The remianing tensor product with the identity
        for j in range(n_particles-i-2):
            temp = np.kron(temp, np.eye(2))

        # Summing all the products
        x_ham.matrix += temp
        
    return x_ham


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

def compute_full_hamiltonian(n_particles, bin_lambda, l_min, dl, kappa):
    """
    It computes the eigenvalues of the Hamiltonian for a given number of particles, a given number of
    lambda values, a given minimum lambda value, a given lambda step and a given number of eigenvalues
    to compute
    
    :param n_particles: number of particles
    :param bin_lambda: number of values of lambda to be considered
    :param l_min: the minimum value of lambda
    :param dl: the step size for the lambda values
    :param kappa: number of eigenvalues to be computed
    :return: The eigenvalues of the Hamiltonian for each value of lambda.
    """
    
    time_spent = []
    
    # Calling the function `pauli` and assigning the output to the variables `sigma_x` and `sigma_z`.
    sigma_x, sigma_z = pauli()
    
    # Creating a matrix of size 2^n_particles \times 2^n_particles$ and assigning the value of 0 to
    # all the elements.
    full_hamiltonian = Matrix(2**n_particles)
    
    # Calculating the eigenvalues of the Hamiltonian for a given value of lambda.
    eigenvalues = []
    for i in range(bin_lambda):
        
        # get the partial time of the computation
        partial = time.time()
            
        # Creating the z-part for the Hamiltonian for the system.
        z_ham = z_hamiltonian(sigma_z, n_particles)
                    
        # Creating the x-part for the Hamiltonian for the system.
        x_ham = x_hamiltonian(sigma_x, n_particles)
                
        # Creating the full Hamitlonian for the Ising model
        full_hamiltonian.matrix = (l_min + (dl*i))*z_ham.matrix + x_ham.matrix
            
        # Finding the eigenvalues and eigenvectors of the full Hamiltonian.
        eigvals, eigvects = np.linalg.eigh(full_hamiltonian.matrix)
        
        # Creating a list of lists, where each lists will contain the eigenvalues.
        temp = []
        for k in range(kappa):
            temp.append(eigvals[k])
        eigenvalues.append(temp)
        
        # Time spent for the current cicles
        time_spent.append(time.time()-partial)
        
        # Printing the progress bar and the remainig time
        progress_bar(i+1, bin_lambda, expected_time(time_spent, i+1, bin_lambda))
        
    return eigenvalues

def plot_results(kappa, n_particles, l_min, l_max, dl, bin_lambda, eigenvalues):
    """
    It plots the evolution of the eigenvalues with the parameter $\lambda$ for a given number of
    particles $ and a given number of eigenvalues $\kappa$
    
    :param kappa: number of eigenvalues to be calculated
    :param n_particles: number of particles in the system
    :param l_min: minimum value of lambda
    :param l_max: the maximum value of lambda
    :param dl: the step size for the parameter $\lambda$
    :param bin_lambda: number of points in the interval [l_min, l_max]
    :param eigenvalues: the eigenvalues of the Hamiltonian for each value of $\lambda$
    """
    
    # Creating a figure with a size of 16x8 inches.
    fig, ax = plt.subplots(figsize=(16,8))    

    # Setting the title of the plot.
    ax.set_title(f'Evolution of ${kappa}$ eigenvalues with $N={n_particles}$ ($\lambda \in [{l_min},{l_max}]$ and $d\lambda$ = {dl})')

    # Setting the label of the x-axis to 'Parameters $\lambda$'
    ax.set_xlabel('Parameters $\lambda$')

    # Setting the label of the y-axis to 'Eigenvalue'.
    ax.set_ylabel('Eigenvalue')

    # Adding a grid to the plot.
    ax.grid()
    
    # Plotting the evolution of the eigenvalues with the parameter $\lambda$ for a given number of
    #     particles $ and a given number of eigenvalues $\kappa$
    for k in range(kappa):

        x_temp = []
        y_temp = []

        for j in range(bin_lambda):

            x_temp.append(l_min + (dl*j))
            y_temp.append(eigenvalues[j][k])

        ax.plot(x_temp, y_temp, label = f"$q.n. = {k}$")

    # A way to make the legend fit in the plot.
    if kappa <=8:

        ax.legend()

    elif 8 < kappa <= 32:

        pos = ax.get_position()

        ax.set_position([pos.x0, pos.y0, pos.width * 0.9, pos.height])

        ax.legend(loc='center right', bbox_to_anchor=(1.25, 0.5))

    # Showing the plot.
    plt.show()

def isign():
    """
    Main function for the computation of the Ising model.
    It computes the full Hamiltonian for a given number of particles
    and a given range of the parameter $\lambda$
    """
    
    # Defining the range of the parameter $\lambda$
    l_min = 0
    l_max = 3

    # Calculating the range of the parameter $\lambda$.
    p_lambda=l_max-l_min 

    # Calculating the number of bins for the parameter $\lambda$.
    bin_lambda = int(abs(p_lambda * 100))
    
    # The number of particles.
    n_particles = 9

    # Calculating the discretization for the parameter $\lambda$.
    dl = p_lambda/(bin_lambda-1)

    # Calculating the number of eigenvalues to extract.
    kappa = 2**n_particles
    
    # Setting the variable `start` to the current time.
    start = time.time()

    # Calling the function `compute_full_hamiltonian` and assigning the output to the variable `eigenvalues`.
    eigenvalues = compute_full_hamiltonian(n_particles, bin_lambda, l_min, dl, kappa)
    
    # Setting the variable `end` to the current time.
    end = time.time()
    
    # Getting the whole time of the computation
    total = end-start
    
    # The above code is calculating and printing the total time it took to run the code.
    print(f"I spent {int(total/60)} minutes and {int(total - int(total/60)*60)} seconds to make the calculation")

    # Plotting the results.
    plot_results(kappa, n_particles, l_min, l_max, dl, bin_lambda, eigenvalues)

if __name__ == '__main__':
    isign()