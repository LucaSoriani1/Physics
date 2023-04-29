# Prime numbers generator
One of the central topics in number theory is prime numbers. Prime numbers are a special set of natural numbers greater than $1$ that are divisible only by $1$ and themselves. The importance of prime numbers resides in various fields such as cryptography, computer science and physics. Generating prime numbers can be a challenging task, especially for large numbers. The study of prime numbers has many applications in fields such as cryptography, coding theory and computer science. For example, the security of many cryptography schemes relies on the factorization of large compound numbers into prime numbers being a difficult problem. The computational complexity for computing prime numbers increases with the order of magnitude of the number $n$. As the value of $n$ increases, the number of iterations required to check whether a number is prime increases, which leads to longer run times (a normal algorithm has complexity $\sim O\left(\sqrt{n}\right)$).

An effective way to check whether a number is prime is to check its divisibility to the square root of $n$, because every factor of a number must be less than or equal to its square root.  If a number $N$ is not prime, then there exists at least one factor, let us call it $F$, other than $1$ and $N$ itself, such that $N = F \cdot Q$.  If $F > \sqrt{N}$, then $Q < \sqrt{N}$, otherwise we would have $F \cdot Q > N$. Therefore, when we check all numbers from $2$ up to $int\left(\sqrt{n}\right) + 1$, if we have not found any factor of $n$ up to that point, then $n$ cannot be divisible by any number greater than its square root.
The distribution of prime numbers is a fundamental problem in number theory, with many important open questions. The prime number theorem, proved by Jacques Hadamard and Charles de la Vallée-Poussin independently in 1896, provides an estimate of the number of prime numbers less than a given number $x$. The theorem states that the number of prime numbers less than $x$ is approximately equal to $\frac{x}{log(x)}$. This means that as $x$ increases, the ratio of primes to all numbers less than $x$ approaches zero.

There are many other important results and conjectures in prime number theory, such as:</p>

* The twin prime number conjecture, which states that there are infinite pairs of prime numbers that differ by 2.
* Goldbach's conjecture, which states that every even integer greater than 2 can be expressed as a sum of two prime numbers.
* The Riemann Hypothesis, a conjecture about the distribution of prime numbers that has many important consequences in number theory and beyond.

The code provides two functions to generate and verify prime numbers, such as:

1. The is_prime() function takes as input an integer n and returns True if n is a prime number and False otherwise. It checks the divisibility of n to the square root of n, since any factor of a number must be less than or equal to its square root.
2. The prime_numbers() function takes as input an integer $x$ and generates the list of prime $x$-numbers using the is_prime() function.

The code uses a simple approach to generate the prime numbers by iterating all numbers starting with $2$ and checking if they are prime using the is_prime() function. Prime numbers are added to the list of primes until the list contains the required number of primes. The main() function generates, by defualt, the first $10000$ primes and calculates the time taken to do so. The list of primes and the execution time are printed in the console.
This will print the first $n=10000$ prime numbers by default and the time taken to execute. Note that the time taken to generate the prime numbers can be significant, especially for large numbers.
The program can be extended and modified for different purposes, such as modeling natural processes or generating random data for testing purposes.

### Functions
* is_prime(n)
This function takes an integer n as input and returns True if n is a prime number, and False otherwise. It checks divisibility of n up to the square root of n, as any factor of a number must be less than or equal to its square root.

* prime_numbers(count)
This function takes an integer count as input and generates a list of the first 'count' prime numbers using the is_prime() function.

### Execution
To run the code, you need to call the main() function, which generates the first 40000 prime numbers and calculates the time taken to do so. The list of primes and the execution time are printed to the console.

Note that the time taken to generate prime numbers can be significant, especially for large numbers. The complexity of generating prime numbers is considered to be O(n * sqrt(n)), which means that the time taken to generate primes increases with the size of n.


### Usage
To use the script, run the following command:

```
python prime.py
```

The script will print the first n=4000 prime numbers by default and the time spent for the execution.

### Modifications and Contributions
This program can be extended and modified to suit different purposes, such as modeling natural processes or generating random data for testing purposes. If you would like to contribute to the program or propose modifications, please send a pull request or contact the author.