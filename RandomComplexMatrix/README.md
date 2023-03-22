# Random Complex Matrix & Adjoint

The adjoint matrix, also known as the Hermitian conjugate or conjugate transpose, of a square matrix is obtained by taking the transpose of the matrix and then taking the complex conjugate of each element. In other words, if $A$ is a matrix, then its adjoint matrix $A^{\dagger}$ is given by $$(A^{\dagger})_{ij} = (A)^{*}_{ji}.$$

In quantum mechanics, the adjoint matrix plays a fundamental role as it is used to describe the complex conjugate of quantum states and operators. In particular, the adjoint matrix of an operator is used to describe the conjugate transpose of the operator, which is used to calculate the expectation value of an observable in quantum mechanics. The adjoint matrix is also used to calculate the probability of measuring a certain value of an observable in a quantum system. Therefore, the adjoint matrix is an important tool for studying the behavior of quantum systems and making predictions about their properties.

This Python program computes the adjoint matrix and trace of a given square matrix. The user can input the values of the matrix manually or generate a random complex matrix. The program saves the results to a text file.

This code is the translation of the Fortran90 code that was the original one.

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

The program will prompt the user to enter the number of rows and columns for the matrix. Then, the user can choose to generate a random complex matrix or input the matrix elements manually.

After the matrix is generated or inputted, the program computes and prints the matrix, its trace, and the adjoint matrix and its trace.

Finally, the program prompts the user whether to save the matrices and their properties to a text file. If the user confirms, the program saves the results to a file named results.txt.

In the Python version, you could insert a not-squared matrix but the trace would not computed.

### Fortran90 code

You will find the original code in Fortran90. To run the .f90 script, you will need a compiler, then some extension in Visual Studio code as:

- C/C++
- Modern Fortran
- Code Runner

### Modifications and Contributions

This program can be extended and modified to suit different purposes, such as modeling natural processes or generating random data for testing purposes. If you would like to contribute to the program or propose modifications, please send a pull request or contact the author.
