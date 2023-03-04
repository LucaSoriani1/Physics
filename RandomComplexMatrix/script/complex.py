import numpy as np  # Import the numpy library

"""
The adjoint matrix, also known as the Hermitian conjugate or conjugate transpose, of a square matrix is obtained by taking the transpose of the matrix and then taking the complex conjugate of each element. In other words, if A is a matrix, then its adjoint matrix A† is given by (A†)ij = (A)ji\*.

In quantum mechanics, the adjoint matrix plays a fundamental role as it is used to describe the complex conjugate of quantum states and operators. In particular, the adjoint matrix of an operator is used to describe the conjugate transpose of the operator, which is used to calculate the expectation value of an observable in quantum mechanics. The adjoint matrix is also used to calculate the probability of measuring a certain value of an observable in a quantum system. Therefore, the adjoint matrix is an important tool for studying the behavior of quantum systems and making predictions about their properties.

This Python program computes the adjoint matrix and trace of a given square matrix. The user can input the values of the matrix manually or generate a random complex matrix. The program saves the results to a text file.

This code is the translation of the Fortran90 code that was the original one.
"""



# The class Matrix has two attributes: dim and mat. The attribute dim is a tuple containing the
# dimensions of the matrix, and the attribute mat is a numpy array filled with zeros.
class Matrix:  # Define a class named Matrix
    def __init__(self, row, col):  # Define the class constructor
        self.dim = (row, col)  # Create an attribute called dim and assign it a tuple containing the dimensions of the matrix
        self.mat = np.zeros((row, col), dtype='complex_')  # Create an attribute called mat and assign it a numpy array filled with zeros
        self.trace = None  # Create an attribute called trace and assign it a None value by default

def create_random_matrix(row, col):  # Define a function to create a random complex matrix
    """
    It creates a random complex matrix with the given dimensions
    
    :param row: The number of rows in the matrix
    :param col: The number of columns in the matrix
    :return: A matrix with random complex numbers
    """
    reals = np.random.rand(row, col)  # Generate a numpy array filled with random numbers between 0 and 1
    imaginary = np.random.rand(row, col)  # Generate another numpy array filled with random numbers between 0 and 1
    matrix = Matrix(row, col)  # Create an instance of the Matrix class with the given dimensions
    for i in range(row):  # Loop through the rows
        for j in range(col):  # Loop through the columns
            matrix.mat[i, j] = reals[i, j] + 1j * imaginary[i, j]  # Assign each element of the matrix a random complex number
    return matrix  # Return the created matrix

def create_user_matrix(row, col):  # Define a function to create a user-defined complex matrix
    """
    It creates a user-defined complex matrix
    
    :param row: The number of rows of the matrix
    :param col: The number of columns in the matrix
    :return: The created matrix
    """
    matrix = Matrix(row, col)  # Create an instance of the Matrix class with the given dimensions
    for i in range(row):  # Loop through the rows
        for j in range(col):  # Loop through the columns
            real_part = float(input(f'Insert the REAL part of the element ({i}, {j}): '))  # Prompt the user to enter the real part of the element
            imaginary_part = float(input(f'Insert the IMAGINARY part of the element ({i}, {j}): '))  # Prompt the user to enter the imaginary part of the element
            matrix.mat[i, j] = real_part + 1j * imaginary_part  # Assign each element of the matrix the user-defined complex number
    return matrix  # Return the created matrix

def calculate_trace(matrix):
    """
    The function takes a `Matrix` object as input and computes the trace of the matrix if it is squared.
    
    :param matrix: The matrix to calculate the trace of
    :return: The trace of the matrix for squared matrix or None for other.
    """
    if matrix.dim[0] != matrix.dim[1]:  # Check if the matrix is not squared
        return None
    trace = np.trace(matrix.mat)  # Compute the trace of the matrix
    matrix.trace = trace  # Store the trace in the `trace` attribute of the `Matrix` object
    return trace

def create_adjoint_matrix(matrix):
    """
    Compute the adjoint of a given matrix and return it as a `Matrix` object.
    The adjoint is also known as the conjugate transpose.
    
    :param matrix: a `Matrix` object
    :return: The adjoint of the given matrix.
    """
    adjoint_matrix = Matrix(matrix.dim[1], matrix.dim[0])  # Create a `Matrix` object with flipped dimensions
    adjoint_matrix.mat = np.conj(matrix.mat.T)  # Compute the adjoint by taking the conjugate transpose
    return adjoint_matrix

def save_to_file(matrix, adjoint_matrix):
    """
    Write the matrices and their properties to a text file.
    The file will contain the original matrix, its dimensions, and its trace (if squared).
    It will also contain the adjoint matrix, its dimensions, and its trace (if squared).
    If the trace cannot be computed, it will be printed as "No trace, not a squared matrix."
    
    :param matrix: a `Matrix` object representing the original matrix
    :param adjoint_matrix: a `Matrix` object representing the adjoint of the original matrix
    """
    filename = 'results.txt'  # Define the filename
    with open(filename, 'w') as f:  # Open the file in write mode
        f.write('The matrix is:\n')
        f.write(str(matrix.mat))  # Write the original matrix
        f.write(f'\nThe dimensions of the matrix is: {matrix.dim[0]} x {matrix.dim[1]}\n')
        trace = calculate_trace(matrix)  # Compute the trace of the original matrix
        if trace is not None:
            f.write(f'The trace of the matrix is: {trace}\n')
        else:
            f.write('No trace, not a squared matrix.\n')

        f.write('\n\nThe ADJOINT Matrix is:\n')
        f.write(str(adjoint_matrix.mat))  # Write the adjoint matrix
        f.write(f'\nThe dimensions of the matrix is: {adjoint_matrix.dim[0]} x {adjoint_matrix.dim[1]}\n')
        trace = calculate_trace(adjoint_matrix)  # Compute the trace of the adjoint matrix
        if trace is not None:
            f.write(f'The trace of the adjoint matrix is: {trace}\n')
        else:
            f.write('No trace, not a squared matrix.\n')
    print(f'The results were saved to {filename}')  # Print a confirmation message

def main():
    """
    It prompts the user to choose whether to create a random or user-defined matrix, then it computes
    the adjoint matrix and its trace and prints them
    :return: the adjoint matrix.
    """
    # Prompting the user to enter the dimensions of the matrix.
    row = int(input('Please, enter the number of rows for the matrix: '))
    col = int(input('Please, enter the number of columns for the matrix: '))

    # Prompting the user to choose whether to create a random or user-defined matrix.
    print('What do you want to do?')
    print('1) RANDOM COMPLEX MATRIX')
    print('2) USER-DEFINED COMPLEX MATRIX')

    choice = int(input())
    if choice == 1:
        matrix = create_random_matrix(row, col)
    elif choice == 2:
        matrix = create_user_matrix(row, col)
    else:
        print('No choice was made')
        return

    # Printing the matrix and its trace.
    print('\nYour Matrix is:\n', matrix.mat)
    trace = calculate_trace(matrix)
    if trace is not None:
        print(f'The trace of the matrix is: {trace}')
    else:
        print('No trace, not a squared matrix')
    print()
    
    # Computing the adjoint of the matrix and its trace and printing them.
    adjoint_matrix = create_adjoint_matrix(matrix)
    print('Your adjoint matrix is:\n')
    print(adjoint_matrix.mat)
    adjoint_trace = calculate_trace(adjoint_matrix)
    if adjoint_trace is not None:
        print(f'The trace of this adjoint matrix is: {adjoint_trace}')
    else: 
        print('No trace, not a square matrix.')
    print()

    # Prompting the user to choose whether to print the matrices to a file.
    confirm = ['y', 'yes', 'si', 'sì']
    choice = input('Do you want to print the matrices to a TXT file? ')
    if choice.lower() in confirm:
        save_to_file(matrix, adjoint_matrix)

if __name__ == '__main__':
    main()