# Fibonacci Sequence and Golden Ratio Calculator
The Fibonacci sequence is a mathematical sequence of numbers in which each number is the sum of the previous two numbers. The sequence begins with $0$ and $1$, and proceeds as follows: $0$, $1$, $2$, $3$, $5$, $8$, $13$, $21$, $34$, $55$, $89$, $144$, and so on. The Fibonacci sequence is found in many natural phenomena, such as the arrangement of leaves on a branch, the spiraling of shells and the branching of trees. The Fibonacci sequence is defined mathematically by the following recursive formula: $$F(n) = F(n-1) + F(n-2),$$ where $F(n)$ is the nth number in the Fibonacci sequence, $F(n-1)$ is the previous number in the sequence, and $F(n-2)$ is the number preceding the previous number.

The Fibonacci sequence is interesting for many reasons. First, it is present in many natural phenomena, such as the arrangement of leaves on a branch, flower petals, mollusk shells, and many others. This is because the succession represents exponential growth in nature. In addition, the Fibonacci sequence has many interesting properties. For example, the ratio of two consecutive numbers in the Fibonacci sequence always approaches the mathematical constant known as the "golden ratio" or "Fibonacci number": $$\phi = \frac{\left(1 + \sqrt{5}\right)}{2} \approx 1,61803398875.$$  This means that if you take one number in the Fibonacci sequence and divide the next, the result always approaches $\phi$. The golden ratio, also known as divine proportion, can be found in nature, art and architecture, and is often associated with beauty and harmony.

Thanks to the value of $\phi$, one can calculate the nth number of the Fibonacci sequence without having to calculate all the previous numbers: $$F(n) = \frac{\left[\phi^n - \left(1-\phi\right)^n\right]}{5}$$

This Python script generates the Fibonacci sequence and calculates the ratio of the last two numbers in the sequence and the difference between this ratio and the golden rational. In the code, the constant 'GOLDEN_RATIO' is defined as described above. 

The function 'fibonacci_sequence(length)' defines the Fibonacci sequence using a for loop. The function takes as argument the length of the sequence you want to generate. A list called 'sequence' containing the first two numbers in the sequence is initialized. A pair of variables $x$ and $y$ is then defined, initialized to $0$ and $1$ respectively, representing the last two numbers in the sequence. Using a for loop, the function generates the rest of the Fibonacci sequence by continuously updating $x$ and $y$ and adding $y$ to the sequence list. Finally, the function returns the sequence list. 


The 'main()' function uses the fibonacci_sequence() function to generate the Fibonacci sequence. The user is prompted to enter the length of the Fibonacci sequence to be generated. Next, the fibonacci_sequence() function is used to generate the Fibonacci sequence and the resulting sequence is printed. The ratio of the last two numbers in the sequence is then calculated and printed. The difference between this ratio and the golden rational, GOLDEN_RATIO, is also calculated and printed.

### Requirements
This script requires Python 3 and the math module, which is included in the standard library.

### Usage
To use the script, run the following command:

```
python fibonacci.py
```

You will be prompted to enter the length of the Fibonacci sequence you wish to calculate. After entering the value, the script will print the first n numbers of the sequence and the ratio between the last two numbers.

Finally, the script will print the value of the golden ratio and the difference between the ratio of the last two Fibonacci numbers and the golden ratio.

### Modifications and Contributions
This program can be extended and modified to suit different purposes, such as modeling natural processes or generating random data for testing purposes. If you would like to contribute to the program or propose modifications, please send a pull request or contact the author.