# Random Trajectories
Random walk, or "random motion" (also known as "drunken man walk"), is a mathematical model that describes the motion of a particle or object moving randomly, without a predefined direction or pattern. In other words, each step or displacement the particle takes is determined randomly, without depending on the previous step. Trajectories characterized by the lack of a clear pattern or direction and can be influenced by various random factors or events. Random trajectories can be observed in various fields, such as physics, biology, finance and computer science, and are often modeled using stochastic processes or probability distributions. Examples of random trajectories are the movements of particles in a gas, the growth patterns of bacterial colonies, fluctuations in stock prices, or the behavior of users on a website. The random walk can be used to describe the behavior of a particle in a liquid or gas, or to analyze stock price movements on the stock market. 

Random motion can be described by a series of steps, or steps, which can be represented by a Markov process. A Markovian process (or Markov chain) is a mathematical model that describes the evolution of a discrete random variable over time, in which the probability of the variable at time $t+1$ depends only on the value taken at time $t$ (and not on all previous values).

A common form of random walk is the "one-dimensional random walk," in which a particle moves along only one dimension (e.g., the $x$-axis). In this case, the position of the particle after n steps can be described by the following formula:$$X_{n} = X_{0} + S_{1} + S_{2} + \ldots + S_{n}$$ where $X_{0}$ is the initial position, $S_{1}$, $S_{2}$, ..., $S_{n}$ are the successive steps that the particle takes randomly, and $X_{n}$ is the position after $n$ steps. 

The random walk is a very simple but extremely versatile model, and it is used in many fields of science and economics. 

*For example, in physics the random walk is often used to describe the motion of particles in a fluid or gas, or to analyze the behavior of polymers.

* In mathematics, random walk has several applications, such as probability theory, game theory, and the theory of stochastic processes. In particular, the random walk is an example of a discrete stochastic process, that is, a process in which the random variable changes only at discrete points in time.

* In finance, random walk is often used to model the behavior of stock prices in the stock market. The random walk theory of prices suggests that stock prices follow a random and unpredictable pattern in the short term, and therefore it is very difficult to make accurate predictions about future price movements.

One of the implications of the random walk is that the final position of the particle is very difficult to predict, since it depends on the outcome of a series of random events. However, it is possible to calculate some statistical properties of the random walk, such as the mean and standard deviation of the final position. 

In this example, the aim was to demonstrate the randomness of trajectories in a Markovian process. The direction of movement and the next position depend only on the current state and not on previous trajectories. This feature makes Markovian processes "memoryless," which means that it is not necessary to keep a record of past positions to determine the next position. Therefore, the trajectory of particles in a Markovian process can take any form. In general, particles can go back and forth between positions, jump to random positions and accumulate random positions.

The project consists of three scripts, all with the same purpose: to test the randomness of trajectories for a Markovian process. Specifically, the scrpits generate a list of random trajectories, where each trajectory is a list of positions, and then plot these trajectories using the Matplotlib library. The programs differ in the dimesion in which the motion occurs, i.e., $1D$, $2D$ or $3D$.
The function 'generate_random_trajectories(NUM_TRAJECTORIES, MAX_STEP)' generates a list of random trajectories. The input parameters 'NUM_TRAJECTORIES' and 'MAX_STEP' specify the number of trajectories to be generated and the maximum number of steps per trajectory, respectively. In the case of the $1D$ script, the function uses a list of possible movements (left or right), creates an empty list to contain all trajectories, and generates 'n_trajs' random trajectories by randomly choosing a direction in which to move at each step. For $2D$, on the other hand, in addition to generating a random displacement of $-1$ or $1$, it also chooses the direction randomly, i.e., whether along the $x$-axis or along the $y$-axis. Similarly, the script for $3D$ works: you choose the random displacement between $-1$ and $1$ and the random direction between $x$, $y$ and $z$.

The 'plot_trajectories(trajectories, MAX_STEP)' takes a list of trajectories and plots them using 'Matplotlib'. The input parameter 'MAX_STEP' specifies the number of steps to be taken for each trajectory. The function plots each trajectory based on their respective positions, using the 'plot()'  function of the Matplotlib library.

Finally, the 'main()' function defines the number of trajectories to be generated and the maximum number of steps for each trajectory, and then calls the generate_random_trajectories() and plot_trajectories() functions to generate and plot the random trajectories.

### Usage
To run the program, open the terminal and type the following command:
```
python random_trajectories.py
```
Five random trajectories will be generated and plotted.

### Dependencies
The program requires the installation of the following Python libraries:

* numpy
* matplotlib
* scipy

The libraries can be installed using pip, for example:
```
pip install -r requirements.txt
```


The default values for NUM_TRAJECTORIES and MAX_STEP are 5 and 90000, respectively. You can change these values in the main function.

### Modifications and Contributions
This program can be extended and modified to suit different purposes, such as modeling natural processes or generating random data for testing purposes. If you would like to contribute to the program or propose modifications, please send a pull request or contact the author.