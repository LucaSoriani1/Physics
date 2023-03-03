import random
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

"""
Random trajectories refer to paths or movements that are unpredictable or irregular. 
They are characterized by a lack of a clear pattern or direction, 
and they can be influenced by various random factors or events. 
Random trajectories can be observed in various fields, such as physics, biology, finance
and computer science, and they are often modeled using stochastic processes or probability distributions.
Examples of random trajectories include the movements of particles in a gas,
the growth patterns of bacterial colonies, the fluctuations in stock prices, or the behavior of users on a website.

In this program, we generate and plot multiple random trajectories. 
Random trajectories are sequences of positions or events that occur randomly or in a random way.
In this program, we generate random trajectories by randomly choosing the direction of movement (left or right) at each step.
We then plot these trajectories using the Matplotlib library.

We then define the list of possible movements and the maximum number of steps for each trajectory. 
Next, we generate n_trajs random trajectories by randomly choosing a direction to move at each step. 
Finally, we plot the trajectories using Matplotlib.

The program provides a simple example of how to generate and plot random trajectories. 
"""


def generate_random_trajectories(NUM_TRAJECTORIES, MAX_STEP):
    """
    > The function generates a list of random trajectories, where each trajectory is a list of positions
    
    :param NUM_TRAJECTORIES: the number of trajectories to generate
    :param MAX_STEP: the maximum number of steps in each trajectory
    :return: A list of lists. Each list is a trajectory.
    """

    # define the list of possible movements: left or right
    movements = [-1, 1]

    # create an empty list to hold all the trajectories
    trajectories = []

    # generate n_trajs random trajectories
    for i in range(NUM_TRAJECTORIES):

        # initialize the position of the trajectory
        x = 0

        # all steps for each trajectories
        trajectory = []

        # for each step, randomly choose a direction to move
        for j in range(MAX_STEP):

            # add the current position to the trajectory
            trajectory.append(x)

            # randomly choose a direction to move
            x += random.choice(movements)

        # add the completed trajectory to the list
        trajectories.append(trajectory)

    return trajectories


def plot_trajectories(trajectories, MAX_STEP):
    """
    It takes a list of trajectories and plots them using Matplotlib
    
    :param trajectories: a list of trajectories, where each trajectory is a list of positions
    :param MAX_STEP: the number of steps to take in each trajectory
    """

    # create a list for the y-axis
    steps = np.arange(MAX_STEP)
    
    # plots the trajectories using Matplotlib.
    fig, ax = plt.subplots()
    ax.set_xlabel('Step')
    ax.set_ylabel('Position')
    ax.set_title('Random Trajectories')
    for traj in trajectories:
        ax.plot(steps, traj)
    plt.show()


def main():
    """
    The function generate_random_trajectories takes in two arguments, NUM_TRAJECTORIES and MAX_STEP,
    and returns a list of trajectories
    """

    # define the number of trajectories to generate
    NUM_TRAJECTORIES = 5

    # define the maximum number of steps for each trajectory
    MAX_STEP = 90000

    # Calling the function generate_random_trajectories and plot_trajectories.
    trajectories = generate_random_trajectories(NUM_TRAJECTORIES, MAX_STEP)
    plot_trajectories(trajectories, MAX_STEP)


if __name__ == '__main__':
    main()