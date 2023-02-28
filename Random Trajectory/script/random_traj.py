import random
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

"""
In this program, we generate and plot multiple random trajectories. 
Random trajectories are sequences of positions or events that occur randomly or in a random way.
In this program, we generate random trajectories by randomly choosing the direction of movement (left or right) at each step.
We then plot these trajectories using the Matplotlib library.

We then define the list of possible movements and the maximum number of steps for each trajectory. 
Next, we generate n_trajs random trajectories by randomly choosing a direction to move at each step. 
Finally, we plot the trajectories using Matplotlib.

The program provides a simple example of how to generate and plot random trajectories. 
It can be extended and modified to suit different purposes, such as modeling natural processes or generating random data for testing purposes.
"""


def generate_random_trajectories(num_trajectories, max_steps):
    """
    This function generates random trajectories by randomly choosing the direction of movement (left or right) at each step.
    It returns a list of all the trajectories generated.
    """

    # define the list of possible movements: left or right
    movements = [-1, 1]

    # create an empty list to hold all the trajectories
    trajectories = []

    # generate n_trajs random trajectories
    for i in range(num_trajectories):

        # initialize the position of the trajectory
        x = 0

        # all steps for each trajectories
        trajectory = []

        # for each step, randomly choose a direction to move
        for j in range(max_steps):

            # add the current position to the trajectory
            trajectory.append(x)

            # randomly choose a direction to move
            x += random.choice(movements)

        # add the completed trajectory to the list
        trajectories.append(trajectory)

    return trajectories


def plot_trajectories(trajectories, max_steps):
    """
    This function plots the trajectories using Matplotlib.
    """

    # create a list for the y-axis
    steps = np.arange(max_steps)
    
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
    This is the main function that runs the program.
    """

    # define the number of trajectories to generate
    num_trajectories = 5

    # define the maximum number of steps for each trajectory
    max_steps = 90000

    # Calling the function generate_random_trajectories and plot_trajectories.
    trajectories = generate_random_trajectories(num_trajectories, max_steps)
    plot_trajectories(trajectories, max_steps)


if __name__ == '__main__':
    main()