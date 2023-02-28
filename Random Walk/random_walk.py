import random
import numpy as np
from matplotlib import pyplot as plt 
from scipy.stats import norm
import math

def simulate_walk(step, limit):
    
    walks = np.zeros((2*step+1, 2))
    
    walks[:, 0] = np.arange(-step, step+1)
    
    for _ in range(limit):
        x = np.random.choice([-1,1], size=step).sum()
        walks[x+step,1] += 1

    walks[:,1] /= limit
    
    return walks

def plot_distribution(step, walks):
    """
    Plot the resulting probability distribution and compare it with a normal distribution.
    """
    mu = 0
    sigma = step**0.5
    
    x1=-step
    x2=step

    z1=((x1-mu)/sigma)-100
    z2=((x2-mu)/sigma)+100

    x1 = np.arange(z1,z2,0.001)
    nd = math.sqrt(2)/math.sqrt(np.pi*sigma*sigma) * np.exp(-(x1*x1)/(2*sigma*sigma))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlim([-step-2, step+2])
    ax.plot(x1, nd, color='red')
    ax.bar(walks[:,0], walks[:,1], 1.8)
    ax.set_xlabel('Position, with l=1')
    ax.set_ylabel('Probability')
    ax.set_title('Probability distribuiton with N = %s step' %step)
    plt.show()
    step=15

STEP = 15
LIMIT = 200000

# Simulazione del random walk e plot della distribuzione di probabilità
walks = simulate_walk(STEP, LIMIT)
np.set_printoptions(precision=5)
print(f"After {LIMIT} cycles, the results are:\n{walks}")
plot_distribution(STEP, walks)