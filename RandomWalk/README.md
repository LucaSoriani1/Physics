# Random Walks

Random walk, or "random motion" (also known as "drunken man walk"), is a mathematical model that describes the motion of a particle or object moving randomly, without a predefined direction or pattern. In other words, each step or displacement the particle takes is determined randomly, without depending on the previous step. Trajectories characterized by the lack of a clear pattern or direction and can be influenced by various random factors or events. Random trajectories can be observed in various fields, such as physics, biology, finance and computer science, and are often modeled using stochastic processes or probability distributions. Examples of random trajectories are the movements of particles in a gas, the growth patterns of bacterial colonies, fluctuations in stock prices, or the behavior of users on a website. The random walk can be used to describe the behavior of a particle in a liquid or gas, or to analyze stock price movements on the stock market. 

Random motion can be described by a series of steps, or steps, which can be represented by a Markov process. A Markovian process (or Markov chain) is a mathematical model that describes the evolution of a discrete random variable over time, in which the probability of the variable at time $t+1$ depends only on the value taken at time $t$ (and not on all previous values).

A common form of random walk is the "one-dimensional random walk," in which a particle moves along only one dimension (e.g., the $x$-axis). In this case, the position of the particle after n steps can be described by the following formula:$$X_{n} = X_{0} + S_{1} + S_{2} + \ldots + S_{n}$$ where $X_{0}$ is the initial position, $S_{1}$, $S_{2}$, ..., $S_{n}$ are the successive steps that the particle takes randomly, and $X_{n}$ is the position after $n$ steps. 

The random walk is a very simple but extremely versatile model, and it is used in many fields of science and economics. 

*For example, in physics the random walk is often used to describe the motion of particles in a fluid or gas, or to analyze the behavior of polymers.

* In mathematics, random walk has several applications, such as probability theory, game theory, and the theory of stochastic processes. In particular, the random walk is an example of a discrete stochastic process, that is, a process in which the random variable changes only at discrete points in time.

* In finance, random walk is often used to model the behavior of stock prices in the stock market. The random walk theory of prices suggests that stock prices follow a random and unpredictable pattern in the short term, and therefore it is very difficult to make accurate predictions about future price movements.

One of the implications of the random walk is that the final position of the particle is very difficult to predict, since it depends on the outcome of a series of random events. However, it is possible to calculate some statistical properties of the random walk, such as the mean and standard deviation of the final position. 

In this example we want to go to demonstrate the distribution of the final position in a one-dimensional random walk. It is assumed that the steps of the random walk process are identically and independently distributed with zero mean and variance $\sigma^{2}$. It is shown that the final position $X_{n}$ follows a Gaussian distribution with zero mean and variance $n\cdot\sigma^{2}$:$$S_{n} \sim N\left(0, n\cdot\sigma^{2}\right)$$

This result can be proved using the central limit theorem, which states that the sum of n independent and identically distributed random variables tends to a Gaussian distribution when n becomes large. This means that the final position of the random walk process will be distributed around zero, with a higher probability of remaining close to the origin than the most extreme positions. The final position distribution is always tends to the same result, regardless of how the Markovian process has evolved over time, making it an equilibrium distribution for the process.


The program calculates the final position distribution of a random walk of a Markovian process. Specifically, the function 'simulate_walk()' simulates a random walk of 'step'> length and repeats it 'limit' times. The final position of each random walk is recorded in the walks array, which has size $(2*step+1, 2)$. The first column of walks is an array of integer values representing the possible end position for the random walk, while the second column keeps track of how often each end position is reached.

After simulating all random walks, the function normalizes the frequencies of the final positions by dividing by the total number of 'limit' simulated.

The 'plot_distribution()' function uses the data obtained from simulate_walk() to plot the distribution of the final position as a bar graph. It also computes the normal distribution for a mean equal to $0$ and a standard deviation equal to $step^{1/2}$. The function also plots the normal distribution as a red line above the bar graph.

Finally, the 'main()' function runs the simulation and plots the resulting distribution for a set number of step steps and a set amount of limit simulations. The program uses the 'numpy' module to perform the vector calculation, the 'matplotlib' module to plot the graphs.

### Usage
To run the program, open the terminal and type the following command:
```
python random_walk.py
```

This will simulate a random walk and plot the resulting probability distribution.

### Dependencies
The program requires the installation of the following Python libraries:

* numpy
* matplotlib
* scipy

The libraries can be installed using pip, for example:
```
pip install -r requirements.txt
```

To use the code, run the main function. The function simulate_walk takes two arguments, step and limit, and returns an array of the positions of the walk and the frequency of each position.

The default values for STEP and LIMIT are 15 and 200000, respectively. You can change these values in the main function.


### Modifications and Contributions
This program can be extended and modified to suit different purposes, such as modeling natural processes or generating random data for testing purposes. If you would like to contribute to the program or propose modifications, please send a pull request or contact the author.