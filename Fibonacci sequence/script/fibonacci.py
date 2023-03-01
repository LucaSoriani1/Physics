import math

"""
This script calculates the Fibonacci sequence up to a given length,
then calculates the ratio between the last two numbers in the sequence
and compares it with the golden ratio.

The Fibonacci sequence is a mathematical sequence of numbers
in which each number is the sum of the two preceding numbers.
The sequence starts with 0 and 1, and goes as follows: 0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, ...

The Fibonacci sequence is found in many natural phenomena,
such as the arrangement of leaves on a stem, the spiral of shells, and the branching of trees.

The golden ratio, also known as the divine proportion,
is a mathematical ratio that is approximately equal to 1.61803398875.
It is found in nature, art, and architecture, and is often associated with beauty and harmony.

The relationship between the Fibonacci sequence and the golden ratio is that the ratio
of any two consecutive numbers in the sequence approaches the golden ratio as the sequence goes to infinity.
In other words, the larger the numbers in the sequence, the closer the ratio between them is to the golden ratio.
"""



# Defining a constant called GOLDEN_RATIO, which is the golden ratio.
GOLDEN_RATIO = (1 + math.sqrt(5)) / 2

def fibonacci_sequence(length):
    """
    We initialize the sequence with the first number in the sequence, then we initialize the first two
    numbers in the sequence, then we generate the sequence up to the desired length by updating the
    current and previous numbers in the sequence and appending the new number to the sequence, and
    finally we return the complete sequence
    
    :param length: The length of the sequence to generate
    :return: The complete sequence
    """

  # Initialize the sequence with the first number in the sequence
    sequence = [0,1]

    # Initialize the first two numbers in the sequence
    x, y = 0, 1

    # Generate the sequence up to the desired length
    for _ in range(length - 2):
        
        # Update the current and previous numbers in the sequence
        x, y = y, x + y
        
        # Append the new number to the sequence
        sequence.append(y)
        
    # Return the complete sequence
    return sequence



def main():
    
    """
    The function takes a length as an argument and returns a list of the first n numbers of the
    Fibonacci sequence
    """
    
   # Prompt the user for the length of the Fibonacci sequence to generate
    length = int(input('How many numbers of the Fibonacci sequence should I calculate? '))

    # Generate the Fibonacci sequence
    sequence = fibonacci_sequence(length)

    # Print the generated sequence
    print(f"The first {length} numbers of the Fibonacci sequence are:\n{sequence}")

    # Calculate the ratio of the last two numbers in the sequence
    x, y = sequence[-1], sequence[-2]
    ratio = x / y

    # Print the ratio of the last two numbers in the sequence
    print(f"The ratio of the last two numbers in the sequence is: {ratio}")

    # Calculate the difference between the calculated ratio and the Golden Ratio
    difference = abs((GOLDEN_RATIO - ratio) / GOLDEN_RATIO) * 100

    # Print the Golden Ratio and the difference between it and the calculated ratio
    print(f"While phi is: {GOLDEN_RATIO}\nThe difference from phi is: {difference}%")



if __name__ == '__main__':
    main()
