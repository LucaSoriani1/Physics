# Random Walks

The random walk is a mathematical model used to describe a system that randomly moves in different directions.The code serves as a tool to visualize and understand the properties of random walks, which have applications in various fields such as finance, physics, and biology. 

This code simulates a random walk by generating a series of random steps, where each step can be either +1 or -1. After each step, the position of the system is accumulated. The resulting probability distribution of the positions is plotted using a histogram and compared to the expected distribution from a normal distribution.

The simulate_walk function takes two arguments, step and limit, and returns an array of the positions of the walk and the frequency of each position.

The plot_distribution function takes the number of steps in the walk and the array of walks as arguments, and plots the probability distribution of the walks as bars and the normal distribution as a line.

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