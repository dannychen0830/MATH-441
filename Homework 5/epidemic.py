# %% imports
import numpy as np
import matplotlib.pyplot as plt
import math

print('imports done!')
# %% Create agents and neighborhood structure (lattice)
# number of people per row 
n = 100 
 # first two entries stores position, third stores state of itself and neighbors
neighborhood = np.zeros([n+2,n+2,5])
neighborhood += 2

#initialize the borders to -1
neighborhood[0,:,:] = -1
neighborhood[n+1,:,:] = -1
neighborhood[:,0,:] = -1
neighborhood[:,n+1,:] = -1

im = plt.pcolormesh(neighborhood[1:n+1,1:n+1:,0],vmin=0,vmax=5)
plt.colorbar(im)
# %% Update Neighbors
# for each cell, index 0 is the agent itself, 1-4 are NWSE neighbors
# state -1 means at the border
def update_neighbors(neighborhood):
    for i in range(1,neighborhood.shape[0]-1):
        for j in range(1,neighborhood.shape[1]-1):
            for k in range(1,5):
                if k == 1:
                    neighborhood[i,j,k] = neighborhood[i,j+1,0]
                if k == 2:
                    neighborhood[i,j,k] = neighborhood[i-1,j,0]
                if k == 3:
                    neighborhood[i,j,k] = neighborhood[i,j-1,0]
                if k == 4:
                    neighborhood[i,j,k] = neighborhood[i+1,j,0]
    return neighborhood

neighborhood = update_neighbors(neighborhood)
print(neighborhood[1,100,:])
# %% Define states and Markov Matrix
# states assignemnts:
# 0: vaccinated, 1: recovered, 2: suseptible, 3, asymptomatic, 
# 4: symptomatic, 5: dead, -1: does not exist 

def markovMatrix(agent):
    #parameters:
    beta = 0.4 # infectious rate
    delta = 0.4 # probability of recovery
    gamma = 0.5 # probability of remaining sick
    eta = 0.2 # probability of asymptomatic
    d = 0.00001 #probability of death when healthy 

    infected_nei = 0
    for i in range(1, 5):
        if agent[i] == 3 or agent[i] == 4:
            infected_nei += 1
    infect_prob = (1-np.exp(-1*beta*infected_nei))

    row0 = [1-d,0,0,0,0,0]
    row1 = [0,1-d,0,delta,delta,0]
    row2 = [0,0,(1-d)*(1-infect_prob),0,0,0]
    row3 = [0,0,(1-d)*(eta)*infect_prob,gamma,0,0]
    row4 = [0,0,(1-d)*(1-eta)*infect_prob,0,gamma,0]
    row5 = [d,d,d,1-gamma-delta,1-gamma-delta,1]

    return np.array([row0,row1,row2,row3,row4,row5])

print(markovMatrix(neighborhood[78,42,:]))
# %% initializing neighborhood
# select 5 infected
xpos = np.random.permutation(n)
ypos = np.random.permutation(n)

for i in range(5):
    print(str(xpos[i]) + ',' + str(ypos[i]))
    neighborhood[xpos[i]+1,ypos[i]+1,0] = 4

neighborhood = update_neighbors(neighborhood)
plt.pcolormesh(neighborhood[1:n+1,1:n+1:,0],vmin=0,vmax=5)
#%% Sample from Markov Matrix
def sampleMM(vec):
    s = np.random.uniform(0,1)
    count = -1
    while s > 0:
        count += 1
        s -= vec[count]
    return count

testvec = np.random.uniform(0,1,5)
testvec /= np.sum(testvec)
print(testvec, sampleMM(testvec))
# %% Main Simulation
time_step = 201
record_time = [0,40,80,120,160,200]

fig,ax = plt.subplots(2,3,constrained_layout=True)

count = 0
for t in range(time_step):
    for i in range(1,n+1):
        for j in range(1,n+1):
            matrix = markovMatrix(neighborhood[i,j,:])
            vector = matrix[:,int(neighborhood[i,j,0])]
            next_state = sampleMM(vector)
            neighborhood[i,j,0] = next_state
    neighborhood = update_neighbors(neighborhood)
    if t == record_time[count]:
        row = math.floor(count/3)
        col = count % 3
        print(str(row) + ',' + str(col))
        ax[row,col].pcolormesh(neighborhood[1:n+1,1:n+1:,0],vmin=0,vmax=5)
        ax[row,col].set_title('t = ' + str(t))
        fig.colorbar(ax[row,col].pcolormesh(neighborhood[1:n+1,1:n+1:,0],vmin=0,vmax=5), ax=ax.flat)
        count += 1
# %%
for i in range(1,n+1):
    for j in range(1,n+1):
        matrix = markovMatrix(neighborhood[i,j,:])
        vector = matrix[:,int(neighborhood[i,j,0])]
        next_state = sampleMM(vector)
        neighborhood[i,j,0] = next_state

plt.pcolormesh(neighborhood[1:n+1,1:n+1:,0],vmin=0,vmax=5)

# %%
