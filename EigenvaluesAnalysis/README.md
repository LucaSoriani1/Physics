# Hermitian and Triangular Eigenvalues analysis

Hermitian matrices are complex square matrices that are symmetric with respect to their self-junction. In other words, if $A$ is a Hermitian matrix, then $A=A^H=(A^T)^{*}$, where $A^H$ represents its self-junction, that is, the conjugate transposition of the matrix $A$.
The eigenvalues of Hermitian matrices are real and are deterministically distributed, which means that once the eigenvalues of the matrix are known, their relative positions are uniquely determined. In particular, the eigenvalues of a Hermitian matrix are uniformly distributed along the real axis. The distribution of the eigenvalues of a Hermitian matrix follows a Poisson distribution:
$$ P(x) = a \cdot x^\alpha \cdot e^{-b\cdot x^\beta}$$
 In other words, if we consider the distance between two consecutive eigenvalues and normalize this distance by the mean of all distances, we obtain a probability distribution that follows Poisson's law. This Poisson distribution is a universal characteristic of chaotic physical systems.
In particular, when considering large Hermitian matrices, their statistical properties are independent of the details of their construction and depend only on the symmetry and size of the matrix itself. 

Triangular matrices are matrices in which all elements above or below the main diagonal are zero; in other words, an upper triangular matrix has all elements below the main diagonal equal to zero, while a lower triangular matrix has all elements above the main diagonal equal to zero.
The eigenvalues of a triangular matrix are the elements on the main diagonal. Specifically, if a matrix is upper triangular, the eigenvalues coincide with the elements on the main diagonal, while if a matrix is lower triangular, the eigenvalues are the same, but in the opposite direction. The eigenvalues of triangular matrices are uniformly distributed along the main diagonal. Regarding the decreasing exponential distribution of the eigenvalues of triangular matrices, this phenomenon is known as the "Wigner phenomenon" and was discovered by Eugene Wigner in 1951. 

In general, random matrices have a spectrum of eigenvalues that follows a well-defined statistical distribution. In the case of triangular matrices, the spectrum follows a decreasing exponential distribution. This means that the probability density of the eigenvalues decreases rapidly as the eigenvalues increase. The Wigner phenomenon has been extensively studied in physics, particularly in relation to the statistical properties of atoms and nuclei. The distribution of eigenvalues is always given by the Poisson distribution where the parameter $\alpha \sim 0$. 

This Python script consists of a class called Matrix that creates a random Hermitian complex matrix of user-specified size, and also creates an upper or lower triangular matrix from the Hermitian matrix. Specifically, the Hermitian complex matrix is generated as a random complex matrix in which the real and imaginary parts of the elements are selected from a uniform distribution in the interval $[-1, 1]$, and then the matrix is made Hermitian by taking the real part of the matrix summed to its conjugate complex transpose divided by $2$.

The code also provides a fit_and_plot_data function that takes as input an array of the values of the eigenvalues of the hermitian matrix and an array of the values of the spaces between the eigenvalues of the upper triangular matrix. The function then plots the distribution of the normalized spaces between the eigenvalues. The function uses the Scipy library to perform the "curve_fit" function to fit a curve to the distribution of spaces between eigenvalues. The fit function is defined as a Poisson distribution, whose specific shape depends on the parameters a, b, c and d passed to the function. Finally, the "fit_and_plot_data" function plots a histogram of the normalized values of the spaces between the eigenvalues and overlays the fit curve. The fit function is specified within the function itself.

Finally, the code executes the function "fit_and_plot_data" on a random Hermitian matrix and the upper or lower triangular matrix derived from the Hermitian matrix. The function is executed 50 times by default. The code plots the results of each run on the same figure to compare the distributions of the spaces between the eigenvalues of the different matrices.

### Dependencies
This code requires the following packages:

* NumPy
* SciPy
* Matplotlib

### Installation

This program requires Python 3 to be installed. All required packages are listed in the requirements.txt file. To install the dependencies, run:

```
pip install -r requirements.txt
```

### Usage

To use this program, navigate to the directory where the matrix.py file is located and run:

```
python matrix.py
```

The program will prompt the user to enter the dimension of the matrix and how many cicles to do.
##### Note: 3000x3000 in 30 cicles spent around 40 minutes to be done.


### Fortran90 code

You will find the original code in Fortran90. To run the .f90 script, you will need a compiler, then some extension in Visual Studio code as:

- C/C++
- Modern Fortran
- Code Runner

Moreover, you need to install the 'lapack' library in order to run it.

### Modifications and Contributions

This program can be extended and modified to suit different purposes, such as modeling natural processes or generating random data for testing purposes. If you would like to contribute to the program or propose modifications, please send a pull request or contact the author.
