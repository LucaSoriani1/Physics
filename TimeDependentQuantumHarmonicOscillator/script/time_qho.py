import time
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import eigh_tridiagonal


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

       
def progress_bar(progress, total, times):
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


def norm(vec, dim, dx):
    """
    It takes a vector, a dimension, and a step size, and returns the vector normalized to 1
    
    :param vec: the vector to be normalized
    :param dim: dimension of the system
    :param dx: the step size in the x direction
    :return: The norm of the vector.
    """
    
    # Calculating the norm of the vector.
    norm_of_vector = 0
    for i in range(len(vec)):
        norm_of_vector = norm_of_vector + (vec[i].real**2 + vec[i].imag**2)*dx
    norm_of_vector = np.sqrt(norm_of_vector)
    
    # Normalizing the vector.
    vec = vec/norm_of_vector
    
    return vec


    
def plot(intspace, dx, T, dt, min, max, bintime, binspace, psi, dk, k_i, time_step, quantum_number):
    """
    It plots the evolution of the wave function in time, and the evolution of the mean position of the
    particle
    
    :param intspace: number of points in the space
    :param dx: the space step
    :param T: total time
    :param dt: time step
    :param min: minimum value of the position
    :param max: the maximum value of the position
    :param bintime: number of time steps
    :param binspace: number of points in the space
    :param psi: the wavefunction
    :param dk: the step in the momentum space
    :param k_i: the k-space grid
    :param time_step: the time step of the simulation
    :param quantum_number: the quantum number of the state you want to plot
    """
    
    time_spent = []
    
    # Setting the variable `start` to the current time.
    start = time.time()
    
    # Just a constant that is added to the wavefunction to get the correct energy.
    offset = quantum_number + 0.5
    
    # Creating a figure with the given parameters.
    fig, ax = plt.subplots(figsize=(16,8))
    ax.set_title('H.O. with time-dependet Hamiltonian evolution (L={}, dx={}, T={}, dt={})'.format(intspace, dx, T, dt))
    ax.set_xlabel('Position')
    ax.set_ylabel('Square modulus of Psi')
    ax.set_xlim([int(min/4), int(max/4)+1])
    ax.set_ylim([0, 1 + offset])
    ax.grid()
    
    # Plotting the wavefunction in the position space.
    ax.plot(np.linspace(min, max, binspace), psi.real**2+psi.imag**2 + offset, label='ground state')
    
    # Calculating the mean position of the particle at the initial time.
    somma = 0
    means=[]
    for i in range(binspace):
        somma = somma + dx*dx*(psi[i].real**2 + psi[i].imag**2)*(i-1)
    means.append(somma)
    
    # Calculating the evolution of the wavefunction in time.
    ii = complex(0,1)
    index = 0
    for i in range(bintime):
        
        partial = time.time()
        
        temp = psi

        potx  = np.array([])
        for j in range(binspace):
            potx = np.append(potx, [((min+(j-1)*dx - i*dt/T)**2)/2])
            temp[j] = np.exp(-ii*potx[j] *dt/2)*temp[j]
        ax.plot(np.linspace(min, max, binspace), potx, label='pot_t={}'.format(i/T*dt))
        
        # Calculating the wavefunction in the momentum space.
        forw = np.fft.fft(temp)
        forw = norm(forw, binspace, dk)
        potk = np.zeros(binspace)
        for j in range(binspace):
            potk[j] = (k_i[j]**2)/2
            forw[j] = np.exp(-ii*potk[j]*dt)*forw[j]
        
        # Calculating the wavefunction in the position space.
        back = np.fft.ifft(forw)
        back = norm(back, binspace, dx)
        for j in range(binspace):
            back[j] = np.exp(-ii*potx[j]*dt/2)*back[j]
        back = norm(back, binspace, dx)
        
        # Assigning the value of back to psi.
        psi = back
        
        # Plotting the wavefunction in the position space.
        ax.plot(np.linspace(min, max, binspace), psi.real**2 + psi.imag**2 + offset, label='q_t={}'.format(i/T*dt))
    
        # Calculating the mean position of the particle at each time step.
        somma = 0
        for j in range(binspace):
            somma = somma + dx*dx*(psi[j].real**2 + psi[j].imag**2)*(j-1)
        means.append(somma)
        
        # Time spent for the current cicles
        time_spent.append(time.time()-partial)
        
        # Printing the progress bar and the remainig time
        progress_bar(i+1, bintime, expected_time(time_spent, i+1, bintime))
        
    # Setting the variable `end` to the current time.
    end = time.time()
    
    # Getting the whole time of the computation
    total = end-start
    # The above code is calculating the total time it took to run the code.
    print(f"I spent {int(total/60)} minutes and {int(total - int(total/60)*60)} seconds to make the calculation")
    
    # show the plot
    plt.show()

    
    def f(x):
        """
        This function takes a number as input and returns the same number as output
        
        :param x: The value to be passed to the function
        :return: The value of x
        """
        return x

    
    # Creating a figure with the given parameters.
    fig, ax = plt.subplots(figsize=(16,8))
    ax.set_title('Evolution of the position')
    ax.set_xlabel('Time step')
    ax.set_ylabel('Mean position')
    ax.grid()

    # Plotting the mean position of the particle and the minimum of the potential.
    ax.plot(time_step, [m + min for m in means], label='Mean position of the particle')
    ax.plot(time_step, f(time_step), label='Minimum of the potential')

    # Plotting the legend of the figure.
    plt.legend()
    
    # Showing the plot.
    plt.show()
    
    
def momentum_space(binspace, eigenvectors, quantum_number, dx, dk):
    """
    It takes the eigenvectors and eigenvalues of the Hamiltonian and returns the momentum space
    wavefunction
    
    :param binspace: number of bins in the space domain
    :param eigenvectors: the eigenvectors of the Hamiltonian matrix
    :param quantum_number: the quantum number of the eigenvector you want to transform
    :param dx: the distance between each point in the position space
    :param dk: the step size in momentum space
    :return: The momentum space wavefunction and the momentum space grid.
    """
    
   # Creating a list from the eigenvector of the Hamiltonian matrix.
    psi = []
    for i in range(binspace):
        psi.append(complex(eigenvectors[i,quantum_number], 0))
        
    # Normalizing the wavefunction.
    psi = norm(psi, binspace, dx)
    
    # Creating a list of the values of the k-momentum.
    k_i = np.zeros(binspace)
    for i in range(1, int((len(k_i)+1)/2)):
        k_i[i] = (i-1)*dk
        k_i[len(k_i)-i] = (i-1)*dk
    k_i[int((len(k_i)+2)/2)] = np.pi/dx
    
    return psi, k_i
    

def tridiagonal_and_eigenproblem(binspace, dx, min, quantum_number):
    """
    The function takes in the number of bins, the bin size, the minimum value of the potential, and the
    quantum number of the eigenvalue you want to find. It then creates a tridiagonal matrix with the
    given parameters, and finds the eigenvalues and eigenvectors of the matrix
    
    :param binspace: The number of bins in the space
    :param dx: The step size of the binspace
    :param min: the minimum value of the x-axis
    :param quantum_number: The quantum number of the eigenvalue you want to find
    :return: The eigenvalues and eigenvectors of the tridiagonal matrix.
    """
    
    # Creating a list of the values of the main diagonal of the tridiagonal matrix.
    main_diagonal = []
    for i in range(binspace):    
        main_diagonal.append((min + (i-1)*dx)**2 + (2/(2*dx)**2))
        
    # Creating a list of the values of the subdiagonal of the tridiagonal matrix.
    subdiagonal = np.full( binspace-1, -1/(2*dx)**2)
        
    # Finding the eigenvalues and eigenvectors of the tridiagonal matrix.
    eigenvalues, eigenvectors = eigh_tridiagonal(main_diagonal, subdiagonal, select='i', select_range=(0,quantum_number + 1), check_finite=False)

    print(f'The eigenvalue is #{quantum_number+1}: {eigenvalues[quantum_number]}')
    
    return eigenvalues, eigenvectors


def time_dependent_harmonic_oscillator():
    """
    We define the minimum and maximum values of the space, and the total space. We define the total
    time. We define the step sizes in the x-direction, the time-direction, and the k-direction. We
    define the quantum number of the state that we want to plot. We define the number of bins in the
    space and in the time. We create a list of time steps. We calculate the eigenvalues and eigenvectors
    of the Hamiltonian. We calculate the wave function in the momentum space. We plot the wave function
    at each time step
    """

    # Defining the minimum and maximum values of the space, and the total space.
    min = -15.0
    max = 15.0
    intspace = max-min
    
    # The total time.
    T = 100
    
    # Dx is the step size in the x-direction, dt is the step size in the time-direction, and dk is the
    # step size in the k-direction. It is called 'discretization'
    dx = 0.05
    dt = 0.05
    dk = np.pi/intspace
    
    # The quantum number of the state that we want to plot.
    quantum_number = 0
    
    # Defining the number of bins in the space and in the time.
    binspace = int(intspace/dx)+1
    bintime = int(T/dt) 
    
    # Creating a list of time steps.
    time_step = []
    for i in range(bintime+1):
        time_step.append((i-1)*dt/T)
    
    # Calculating the eigenvalues and eigenvectors of the tridiagonal matrix.
    eigenvalues, eigenvectors = tridiagonal_and_eigenproblem(binspace, dx, min, quantum_number)
    
    # Calculating the wave function in the momentum space.
    psi, k_i = momentum_space(binspace, eigenvectors, quantum_number, dx, dk)
    
    # Plotting the wave function at each time step.
    plot(intspace, dx, T, dt, min, max, bintime, binspace, psi, dk, k_i, time_step, quantum_number)
    
if __name__ == '__main__':
    time_dependent_harmonic_oscillator()

