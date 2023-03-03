import random
import numpy as np
from matplotlib import pyplot as plt 
from scipy.stats import norm
import math

"""
The random walk is a mathematical model used to describe a system that randomly moves in different directions. 
The provided code simulates a random walk by generating a series of random steps, 
each of which can be either +1 or -1, and accumulating the position of the system after each step.
The code then plots the resulting probability distribution of the positions after a certain number of steps,
using a histogram and comparing it to the expected distribution from a normal distribution.
The code serves as a tool to visualize and understand the properties of random walks,
which have applications in various fields such as finance, physics, and biology.
"""


def simulate_walk(step, limit):
    """
    We simulate a random walk by choosing -1 or 1 with equal probability 'step' times, then summing up
    the results to get the final position of the walk
    
    :param step: the number of steps in the random walk
    :param limit: the number of simulations to run
    :return: An array of the positions of the walk and the frequency of each position.
    """
    # Initialize the array of walks with zeros, with dimensions (2*step+1, 2)
    walks = np.zeros((2*step+1, 2))
    
    # Set the first column of the array to be a range of values from -step to step+1
    walks[:, 0] = np.arange(-step, step+1)
    
    # Simulate a random walk by choosing -1 or 1 with equal probability 'step' times, 
    # then summing up the results to get the final position of the walk.
    # Update the corresponding element of the walks array by adding 1 to the second column (frequency)
    for _ in range(limit):
        x = np.random.choice([-1,1], size=step).sum()
        walks[x+step,1] += 1

    # Normalize the frequency of each position by dividing by the total number of simulations 'limit'
    walks[:,1] /= limit
    
    return walks

def plot_distribution(step, walks):
    """
    It takes a step number and a list of walks, and plots the probability distribution of the walks as
    bars, and the normal distribution as a line
    
    :param step: The number of steps in the walk
    :param walks: The array of walks
    """
    
    # Calculate the mean and standard deviation of the walk
    mu = 0
    sigma = step**0.5
    
    # Set the range of x values for the normal distribution plot
    x1=-step
    x2=step

    # Set the range of z values for the normal distribution plot
    # The formula for z is (x - mu) / sigma, where mu is the mean and sigma is the standard deviation
    z1=((x1-mu)/sigma)-100
    z2=((x2-mu)/sigma)+100

    # Set the range of x values to be plotted for the normal distribution
    x1 = np.arange(z1,z2,0.001)

    # Calculate the normal distribution values for each x value
    nd = math.sqrt(2)/math.sqrt(np.pi*sigma*sigma) * np.exp(-(x1*x1)/(2*sigma*sigma))

    # Create a new plot
    fig = plt.figure()

    # Add a subplot to the plot
    ax = fig.add_subplot(111)

    # Set the x-axis limits
    ax.set_xlim([-step-2, step+2])

    # Plot the normal distribution in red
    ax.plot(x1, nd, color='red')

    # Plot the simulated walk distribution as bars
    ax.bar(walks[:,0], walks[:,1], 1.8)

    # Add labels to the x and y axes
    ax.set_xlabel('Position, with l=1')
    ax.set_ylabel('Probability')

    # Add a title to the plot
    ax.set_title('Probability distribuiton with N = %s step' %step)

    # Display the plot
    plt.show()

def main():
    """
    The function `simulate_walk` takes two arguments, `STEP` and `LIMIT`, and returns a list of lists,
    where each list is a random walk of length `STEP`
    """
    
    # Set the number of steps and the limit of simulations
    STEP = 15
    LIMIT = 200000

    # Simulate the random walk and plot the resulting distribution
    walks = simulate_walk(STEP, LIMIT)
    np.set_printoptions(precision=5)
    plot_distribution(STEP, walks)
    
if __name__ == '__main__':
    main()