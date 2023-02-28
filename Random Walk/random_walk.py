import random
import numpy as np
from matplotlib import pyplot as plt 
from scipy.stats import norm
import math




step=int(input('How many step do you want?   '))

ris = np.zeros((2*step+1,2))

jump=(-1,1)


incr = 0
for x in range (2*step+1):
    ris[x][0] = -step + incr
    incr += 1

limit=200000
N=1
while N<=limit:
    x=0
    for y in range (step):
        random_step = random.choice(jump)
        x += random_step
        
    ris[x+step][1] += 1
    
    N+=1

ris[:,1] = ris[:,1]/limit

np.set_printoptions(precision=5)

print('After', limit, 'cycles, the results are: \n', ris)

prob = 0
for g in range (2*step+1):
     prob += ris[g][1]

print('') 
print('The total probability is:', prob)

mu = 0
sigma = math.sqrt(step)
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
ax.bar(ris[:,0], ris[:,1], 1.8)
ax.set_xlabel('Position, with l=1')
ax.set_ylabel('Probability')
ax.set_title('Probability distribuiton with N = %s step' %step)
plt.show()





