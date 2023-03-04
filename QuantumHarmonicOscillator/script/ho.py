import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh_tridiagonal
from tabulate import tabulate

def print_results(eigvals):
    """
    It takes a list of eigenvalues and prints a table of the expected eigenvalues, the computed
    eigenvalues, and the percentage difference between the two
    
    :param eigvals: the eigenvalues of the matrix
    """
    table = []
    
    # Creating a table of the expected eigenvalues, the computed eigenvalues, and the percentage
    # difference between the two.
    for i in range(len(eigvals)):
        expected = round(i+0.5, 3)
        computed = round(eigvals[i], 3)
        
        perc = str(round((expected/computed) * 100, 2)) + '%'
        table.append([i, expected, computed, perc])
        
    # It prints a table of the expected eigenvalues, the computed eigenvalues, and the percentage
    # difference between the two.
    print (tabulate(table, headers=["n", "Expected", "Computed", "Percentage"]))
    

def quantum_harmonic_oscillator(x_min, x_max, discretization_dx, n_of_eigvs):
    """
    It solves the eigenvalue problem for the quantum harmonic oscillator in 1D
    
    :param x_min: the minimum value of the x-axis
    :param x_max: the maximum value of x
    :param discretization_dx: the step size of the discretization of the interval
    :param n_of_eigvs: number of eigenvalues to be calculated
    """
    
    # Calculating the interval of the x-axis.
    interv = x_max - x_min
    
    # The mean of the interval of the x-axis.
    mean = (x_max + x_min) / 2
    
    # Calculating the number of steps in the interval.
    steps = int(interv / discretization_dx) + 1
    
    # Creating a matrix of the potential energy of the quantum harmonic oscillator.
    pot = np.zeros(steps)
    diag = np.zeros(steps)    
    for i in range(steps):

        pot[i] = ((x_min - mean) + (i - 1) * discretization_dx) ** 2

        diag[i] = pot[i] + 2 / (2 * discretization_dx) ** 2
    
    subdiag = np.zeros(steps - 1)
    for i in range(steps - 1):
        subdiag[i] = -1 / (2 * discretization_dx) ** 2
    
    # Solving the eigenvalue problem for the quantum harmonic oscillator in 1D.
    eigv = np.zeros(steps)
    yi = np.zeros((steps, n_of_eigvs))
    work = np.zeros(5 * steps)
    iwork = np.zeros(5 * steps, dtype=np.int32)
    ifail = np.zeros(steps, dtype=np.int32)
    
    eigv, yi = eigh_tridiagonal(diag, subdiag, select='i', select_range=(0, n_of_eigvs-1), check_finite=False)
    
    # It prints a table of the expected eigenvalues, the computed eigenvalues, and the percentage
    # difference between the two.
    print_results(eigv)
        
    # Plotting the eigenvalues and eigenvectors of the quantum harmonic oscillator.
    
    x_minimo = 0.0
    
    # Finding the minimum value of the potential energy of the quantum harmonic oscillator.
    for i in range(steps):
        if pot[i] < x_minimo:
            x_minimo = pot[i]
            
    # Calculating the difference between the maximum and minimum eigenvalues.
    de = (eigv[n_of_eigvs - 1] - eigv[0]) / float(n_of_eigvs - 1)
    
    # Plotting the eigenvalues and eigenfunctions of the quantum harmonic oscillator.
    plt.figure(figsize=(16,8))
    plt.title('Eigenfunction of harmonic oscillator in 1D')
    plt.ylabel('Energies')
    plt.xlabel('Position')
    plt.ylim(x_minimo, eigv[n_of_eigvs - 1] + de)
    plt.xlim(int(x_min), int(x_max))
    de = 1 / np.sqrt(discretization_dx)
    for i in range(n_of_eigvs):
        plt.plot(x_min + np.arange(steps) * discretization_dx, eigv[i] + de * yi[:, i])
    plt.plot(x_min + np.arange(steps) * discretization_dx, pot)
    plt.show()

def main():
    """
    It takes in the minimum and maximum values of the x-axis, the discretization step size, and the
    number of eigenvalues to be calculated, and then it calculates the eigenvalues and eigenvectors of
    the quantum harmonic oscillator
    """
    
    # Defining the interval of the x-axis.
    x_min = -6
    x_max = 6
    
    # The step size of the discretization of the interval.
    discretization_dx = 0.0001
    
    # The number of eigenvalues to be calculated.
    n_of_eigvs=10
    
    # Solving the eigenvalue problem for the quantum harmonic oscillator in 1D.
    quantum_harmonic_oscillator(x_min, x_max, discretization_dx, n_of_eigvs)
    
if __name__ == "__main__":
    main()