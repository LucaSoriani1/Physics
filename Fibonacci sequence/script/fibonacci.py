import math

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
