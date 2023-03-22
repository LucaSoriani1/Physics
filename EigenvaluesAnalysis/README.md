# Hermitian and Triangular Eigenvalues analysis

Hermitian matrices are complex square matrices that are symmetric with respect to their self-adjoint. In other words,
if A is a Hermitian matrix, then A=A^H, where A^H represents its self-adjoint, which is the conjugate transpose of the matrix A.
The eigenvalues of Hermitian matrices are real and are deterministically distributed, which means that once the eigenvalues of the matrix are known, their relative position is uniquely determined. In particular, the eigenvalues of a Hermitian matrix are uniformly distributed along the real axis. The distribution of eigenvalues of a Hermitian matrix follows a Poisson distribution. In other words, if the distance between two consecutive eigenvalues is considered and this distance is normalized by the average of all distances, a probability distribution that follows the Poisson law is obtained. This Poisson distribution is a universal characteristic of disordered physical systems.
In particular, if large Hermitian matrices are considered, their statistical properties are independent of the details of their construction
and depend only on the symmetry and dimension of the matrix itself. 

Triangular matrices are matrices in which all elements above or below the main diagonal are zero, meaning an upper triangular matrix has all elements below the main diagonal equal to zero, while a lower triangular matrix has all elements above the main diagonal equal to zero.
The eigenvalues of a triangular matrix are the elements on the main diagonal. In particular, if a matrix is upper triangular, then the eigenvalues coincide with the elements on the main diagonal, while if a matrix is lower triangular, the eigenvalues are the same, but in the opposite direction. The eigenvalues of triangular matrices are uniformly distributed along the main diagonal. As for the decreasing exponential distribution of eigenvalues of triangular matrices, this phenomenon is known as the "Wigner phenomenon" and was discovered by Eugene Wigner in 1951. 
In general, random matrices have an eigenvalue spectrum that follows a well-defined statistical distribution. In the case of triangular matrices,
the spectrum follows a decreasing exponential distribution. This means that the probability density of eigenvalues decreases rapidly as the eigenvalues themselves increase. The Wigner phenomenon has been widely studied in physics, particularly in relation to the statistical properties of atoms and nuclei.

This code is a Python script that analyzes random matrices, and is designed to generate two triangular matrices from a random Hermitian matrix. The distribution of the eigenvalues of the Hermitian matrix is plotted, and the normalized spacings of the eigenvalues are shown in a histogram. The function curve_fit is used to fit the histogram to a general distribution with exponential shape.

This code is the translation of the Fortran90 code that was the original one.

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
