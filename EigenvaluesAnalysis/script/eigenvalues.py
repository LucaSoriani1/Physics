import numpy as np
import scipy.linalg as la
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import time

"""
Hermitian matrices are complex square matrices that are symmetric
with respect to their self-adjoint. In other words,
if A is a Hermitian matrix, then A=A^H, where A^H represents its self-adjoint,
which is the conjugate transpose of the matrix A.
The eigenvalues of Hermitian matrices are real and are deterministically distributed,
which means that once the eigenvalues of the matrix are known,
their relative position is uniquely determined.
In particular, the eigenvalues of a Hermitian matrix are uniformly distributed along the real axis.
The distribution of eigenvalues of a Hermitian matrix follows a Poisson distribution.
In other words, if the distance between two consecutive eigenvalues is considered
and this distance is normalized by the average of all distances,
a probability distribution that follows the Poisson law is obtained.
This Poisson distribution is a universal characteristic of disordered physical systems.
In particular, if large Hermitian matrices are considered,
their statistical properties are independent of the details of their construction
and depend only on the symmetry and dimension of the matrix itself.

Triangular matrices are matrices in which all elements above or below the main diagonal are zero,
meaning an upper triangular matrix has all elements below the main diagonal equal to zero,
while a lower triangular matrix has all elements above the main diagonal equal to zero.
The eigenvalues of a triangular matrix are the elements on the main diagonal.
In particular, if a matrix is upper triangular,
then the eigenvalues coincide with the elements on the main diagonal,
while if a matrix is lower triangular, the eigenvalues are the same,
but in the opposite direction.
The eigenvalues of triangular matrices are uniformly distributed along the main diagonal.
As for the decreasing exponential distribution of eigenvalues of triangular matrices,
this phenomenon is known as the "Wigner phenomenon" and was discovered by Eugene Wigner in 1951.
In general, random matrices have an eigenvalue spectrum
that follows a well-defined statistical distribution. In the case of triangular matrices,
the spectrum follows a decreasing exponential distribution.
This means that the probability density of eigenvalues decreases rapidly
as the eigenvalues themselves increase. The Wigner phenomenon has been widely studied in physics,
particularly in relation to the statistical properties of atoms and nuclei.
"""

# The class Matrix creates a random hermitian matrix, and then creates two triangular matrices from it
class Matrix():
    
    def __init__(self, dim):
        """
        The function takes a dimension as an input and creates a random hermitian matrix of that
        dimension, a upper triangular and a lower triangualr from the hermitian
        
        :param dim: the dimension of the matrix
        """
        
        # Assigning the value of the parameter `dim` to the attribute `dim` of the class.
        self.dim = dim
        
        # Creating a random complex matrix.
        re = 2*np.random.rand(self.dim, self.dim) - np.ones((self.dim, self.dim))
        im = 2*np.random.rand(self.dim, self.dim) - np.ones((self.dim, self.dim))
        complex_matrix = re + 1j*im
        
        # Creating a hermitian matrix from a complex matrix.
        self.hermitian = 0.5*(complex_matrix + complex_matrix.T.conj())
        
        # Creating a upper and lower triangular matrix from the hermitian matrix.
        self.triangularU = np.triu(self.hermitian)
        self.triangularL = np.tril(self.hermitian)



def fit_and_plot_data(h_array, d_array, title, cicles,  bin=250, x_range=[-0.05,5]):
    """
    It takes the eigenvalues of a Hermitian and a Triangular matrix, and it plots the distribution of
    the normalized spacings of the eigenvalues
    
    :param h_array: the array of eigenvalues of the hermitian matrix
    :param d_array: the array of spacings between eigenvalues of the upper triangular matrix
    :param title: the title of the plot
    :param cicles: number of iterations
    :param bin: number of bins, defaults to 250 (optional)
    :param x_range: the range of the x-axis
    """

    # The number of iterations that the `curve_fit` function will do.
    iteration = len(h_array)
    
    def f(x, a, b, c, d):
        """
        A function that takes in 4 parameters and returns a value.
        It define a general distribution with exponential shape. 
        
        :param x: the x-axis data
        :param a: the amplitude of the curve
        :param b: the exponent of the power law
        :param c: the number of clusters
        :param d: the number of data points
        :return: the value of the function f(x)
        """
        return a*(x**b)*np.exp(-c*(x**d))
    
    # It creates a figure with two subplots, one on the left and one on the right.
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 6))

    # HERMITIAN MATRIX PLOT
    
    # Creating a list of the absolute values of the eigenvalues of the hermitian matrix.
    h_values = [abs(float(x)) for x in h_array if float(x) >= 0]
    
    # Finding the maximum and minimum values of the eigenvalues of the hermitian matrix.
    h_maxi = max(h_values)
    h_mini = min(h_values)

    # Calculating the width of each bin.
    h_dimensione = (h_maxi - h_mini) / bin
    
    # It creates an array of zeros of length `bin`.
    h_val = np.zeros(bin)

    # Counting the number of eigenvalues in each bin.
    for i in range(bin):
        for j in range(len(h_values)):
            if h_values[j] >= (h_mini + i * h_dimensione) and h_values[j] < h_mini + (i+1) * h_dimensione:
                h_val[i] += 1

    # The sum of the values of the array `h_val`.
    h_b = np.sum(h_val)
    
    # Creating a list of values that are the midpoints of the bins.
    h_x = [h_mini + (float(i) + 0.5) * h_dimensione for i in range(len(h_val))]
    
    # Calculating the normalized frequencies of the eigenvalues of the hermitian matrix.
    h_y = [h_val[i] / (h_b * h_dimensione) for i in range(len(h_val))]
    
    # A loop that tries to fit the data with the function `f` until it finds a good fit.
    i=1
    while True:
        try:
            print(f"Tentative #{i}: try to fit hermitian with maxfev = {iteration}")
            h_popt, h_pcov = curve_fit(f, xdata=h_x, ydata=h_y, maxfev = iteration)
            break
        except:
            iteration += len(h_array)
            i+=1
            
    # Assigning the values of the array `h_popt` to the variables `h_a`, `h_b`, `h_c`, `h_d`.
    h_a, h_b, h_c, h_d = h_popt
    
    # Calculating the standard deviation of the parameters of the fit.
    a_err_h, b_err_h, c_err_h, d_err_h = np.sqrt(np.abs(np.nan_to_num(np.diag(h_pcov))))
    
    
    # Generate the plot for the hermitian matrix.
    ax1.bar(h_x, h_y, width=h_dimensione, color='g', alpha=0.5, label='Frequencies normalized spacings')
    ax1.plot(h_x, f(h_x, h_a, h_b, h_c, h_d), 'b-', label='Fit of Distribution')
    ax1.set_xlabel('Values')
    ax1.set_ylabel('Frequencies')
    ax1.set_title('Fit of distribution of normalized spacings of {}x{} Hermitian random matrix'.format(title, title))
    ax1.legend(loc='upper right')
    ax1.set_xlim(x_range)
    
    # Triangular MATRIX PLOT
    
    # Creating a list of the absolute values of the eigenvalues of the upper triangular matrix.
    d_values = [abs(float(x)) for x in d_array if float(x) >= 0]

    # Finding the maximum and minimum values of the list d_values.
    d_maxi = max(d_values)
    d_mini = min(d_values)

    # Calculating the width of each bin.
    d_dimensione = (d_maxi - d_mini) / bin
    
    # It creates an array of zeros of length `bin`. 
    d_val = np.zeros(bin)

    # Counting the number of eigenvalues in each bin.
    for i in range(bin):
        for j in range(len(d_values)):
            if d_values[j] >= (d_mini + i * d_dimensione) and d_values[j] < d_mini + (i+1) * d_dimensione:
                d_val[i] += 1
    # The sum of the values of the array `d_val`.
    d_b = np.sum(d_val)
    
    # Creating a list of values that are the midpoints of the bins.
    d_x = [d_mini + (float(i) + 0.5) * d_dimensione for i in range(len(d_val))]
    
    
    # Calculating the normalized frequencies of the eigenvalues of the hermitian matrix.
    d_y = [d_val[i] / (d_b * d_dimensione) for i in range(len(d_val))]
    
    iteration = len(d_array)
    
    # A loop that tries to fit the data with the function `f` until it finds a good fit.
    i=1
    while True:
        try:
            print(f"Tentative #{i}: try to fit triangular with maxfev = {iteration}")
            d_popt, d_pcov = curve_fit(f, xdata=d_x, ydata=d_y, p0=[10, 0.5, 0.5, 0.5], maxfev=iteration)
            break
        except:
            iteration += len(h_array)
            i+=1
    # Assigning the values of the array `h_popt` to the variables `d_a`, `d_b`, `d_c`, `d_d`.
    d_a, d_b, d_c, d_d = d_popt
    
    # Calculating the standard deviation of the parameters of the fit.
    a_err_d, b_err_d, c_err_d, d_err_d = np.sqrt(np.diag(d_pcov))

    # Generate the plot for the hermitian matrix.
    ax2.bar(d_x, d_y, width=d_dimensione, color='g', alpha=0.5, label='Frequencies normalized spacings')
    ax2.plot(d_x, f(d_x, d_a, d_b, d_c, d_d), 'b-', label='Fit of Distribution')
    ax2.set_xlabel('Values')
    ax2.set_ylabel('Frequencies')
    ax2.set_title('Fit of distribution of normalized spacings of {}x{} Upper Triangular random matrix'.format(title, title))
    ax2.legend(loc='upper right')
    ax2.set_xlim(x_range)

    # Plot the graphs
    plt.show()

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

def main():
    """
    It takes the dimension of the matrices and the number of cycles to do, then it creates the matrices,
    calculates the eigenvalues and the spacings, then it plots the spacings and fits them with the
    given distribution
    """
    # Asking the user to input the dimension of the matrix and the number of cycles.
    dim = int(input('What is the values of the dimension N? '))
    cicles = int(input('How many cycles do you want to do? '))
    
    # Defining the empy array used for the calculation
    h_means = []
    t_means = []
    times = []
    
    # Taking the starting time
    time2 = time.time()
    
    # Generating an hermitian and upper triangular matrices 
    # and then calculating the eigenvalues of theese.
    for l in range(cicles):
        
        # Creating a variable called time1 and assigning it the value of the current time.
        # Starting time of the current cicle
        time1 = time.time()
        
        # Creating an object Matrix of size dim x dim.
        matrix = Matrix(dim)
        
        # Finding the eigenvalues of the hermitian matrix.
        eigvals_hermitian, _ = np.linalg.eigh(matrix.hermitian)
        
        # The above code is calculating the difference between the eigenvalues of the hermitian matrix and
        # then dividing it by the mean of the difference in order to see the distribution.
        difference_eigenvals_hermitian = np.diff(eigvals_hermitian) 
        mean_hermitian = np.mean(difference_eigenvals_hermitian)        
        h_means = np.hstack((h_means, difference_eigenvals_hermitian/mean_hermitian))
        
            
        # Finding the eigenvalues of the triangular matrix.
        eigval_triangular, _ = np.linalg.eigh(matrix.triangularU)
        
        # The above code is calculating the difference between the eigenvalues of the triangular matrix and
        # then dividing it by the mean of the difference in order to see the distribution.
        difference_eigenvals_triangular = np.diff(eigval_triangular) 
        mean_triangular = np.mean(difference_eigenvals_triangular)        
        t_means = np.hstack((t_means, difference_eigenvals_triangular/mean_triangular))
        
        # Taking the end time for the current cicle
        end = time.time()
        
        # Compute the time for the current cicle
        time1 = end - time1
        
        # Populating the times array in order to estimate the remaining time
        times.append(time1)
        
        # Printing a progress bar.
        progress_bar(l+1, cicles, expected_time(times, l+1, cicles))
        
    # Taking the final time of the whole computation
    end = time.time()
    
    # Calculating the total time of the whole computation
    time2 = end - time2
    
    # Printing the time it takes to run the code.
    print(f"TOTAL TIME: {int(time2/60)} minutes and {int((time2-int(time2/60))*60)} seconds")
    
    # Fitting and plotting the results
    fit_and_plot_data(h_means, t_means, title=dim, cicles=cicles)

if __name__ == '__main__':
    main()