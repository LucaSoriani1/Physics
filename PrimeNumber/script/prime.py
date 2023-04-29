import time 

"""
Prime numbers are an essential concept in mathematics,
and they are natural numbers greater than 1 that are divisible only by 1 and themselves.
Prime numbers play a crucial role in number theory and have applications in various fields like cryptography,
computer science, and physics. However, finding prime numbers can be a challenging task, especially for large numbers.
The complexity of generating prime numbers is considered to be O(n * sqrt(n)),
hich means that the time taken to generate primes increases with the size of n.
This is due to the fact that for a given number n, we only need to check its divisibility
up to the square root of n because any factor of a number must be less than or equal to its square root.
As the value of n increases, the number of iterations required to check if a number is prime increases,
leading to longer execution times.

This code defines two functions to generate and check prime numbers.
The first function, is_prime(n), takes an integer n as input and
returns True if n is a prime number, and False otherwise.
It checks divisibility of n up to the square root of n,
as any factor of a number must be less than or equal to its square root.
The second function, prime_numbers(count), takes an integer count as input and
generates a list of the first 'count' prime numbers using the is_prime() function.

The main function generates the first 40000 prime numbers and calculates the time taken to do so.
The code uses a simple approach to generate prime numbers by iterating through all numbers from 2 and
checking if they are prime using the is_prime() function.
Prime numbers are added to the list of primes until the list contains the required number of primes.
The time taken to generate prime numbers can be significant, especially for large numbers.
The complexity of generating prime numbers is considered to be O(n * sqrt(n)),
which means that the time taken to generate primes increases with the size of n.
"""


#Function to check if a number is prime
def is_prime(n):
    """
    This function checks if a number is prime or not.
    
    The reason we only need to check up to the square root of n is that factors
    of a number must be less than or equal to its square root.
    If a number N is not prime, then there exists at least one factor, 
    call it F, that is different from 1 and N itself such that N = F * Q. 
    If F > sqrt(N), then Q < sqrt(N), otherwise we would have F * Q > N.
    Therefore, when we check all numbers from 2 up to int(n ** 0.5) + 1, 
    if we have not found any factors of n up to that point, 
    then n cannot be divisible by any number greater than its square root.

    Args:
    - n (int): the number to be checked

    Returns:
    - bool: True if the number is prime, False otherwise
    """
    # A number less than 2 is not prime
    if n < 2:
        return False
    
    # Check divisibility up to the square root of the number
    for i in range(2, int(n ** 0.5) + 1):
        if n % i == 0:
            return False
    # If the number is not divisible by any number up to its square root, it is prime
    return True

#Function to compute the first n prime numbers
def prime_numbers(count):
    """
    This function generates a list of prime numbers of a given count.

    sql
    Copy code
    Args:
    - count (int): the number of prime numbers to generate

    Returns:
    - list: a list of the first 'count' prime numbers
    """
    # Initialize the list of primes with the first prime number (1)
    primes = [1]
    # Start checking numbers from 2 and keep adding to the list until it has the required number of primes
    n = 2
    while len(primes) < count:
        if is_prime(n):
            primes.append(n)
        n += 1
    # Return the list of primes
    return primes

#Main function
def main():
    # Generate the first 10000 prime numbers
    count = 10000
    start = time.time()

    primes = prime_numbers(count)
    
    end = time.time() 
    execution_time = end-start       
    # Print the list of primes
    print(primes)
    print("Execution time:", execution_time, 's')

if __name__ == '__main__':
    main()
