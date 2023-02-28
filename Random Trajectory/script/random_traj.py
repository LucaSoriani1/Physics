import random
import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import norm
import math

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

# define the list of possible movements: left or right
jump = (-1, 1)

# define the maximum number of steps for each trajectory
max_step = 90000

# create a list for the y-axis
steps = list(range(max_step))

# define the number of trajectories to generate
n_trajs = 5

# create an empty list to hold all the trajectories
trajs = []

# generate n_trajs random trajectories
for i in range(n_trajs):
    # initialize the position of the trajectory
    x = 0
    temp = []
    # for each step, randomly choose a direction to move
    for j in range(max_step):
        # add the current position to the trajectory
        temp.append(x)
        
        # randomly choose a direction to move
        x += random.choice(jump)
        
    # add the completed trajectory to the list
    trajs.append(temp)

# create a new figure and axis for plotting
fig = plt.figure()
ax = fig.add_subplot(111)
# set the x and y axis labels and title
ax.set_xlabel('Step')
ax.set_ylabel('Position')
ax.set_title('Random trajectories')
# plot each trajectory
for el in trajs:
    ax.plot(steps, el)

# display the plot
plt.show()