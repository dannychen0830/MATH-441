#%% 
import matplotlib.pyplot as plt
import numpy as np 
import matplotlib.colors
import math

omega = [120,120] # dimension of the forest field 
d = 0.5 # density of the forest
n = 200 # number of unit cells
#%% 
N = int(d*n*n) # number of trees 

trees = np.zeros([N,2]) # stores the position of each tree
for i in range(N):
    # uniformly place each tree
    trees[i,:] = [np.random.uniform(0,omega[0]), np.random.uniform(0,omega[1])]

plt.scatter(trees[:,0], trees[:,1], c = 'g', s = 1)
plt.xlim(0, omega[0])
plt.ylim(0, omega[1])
#%% 
# keep track of the states of the trees
unburnt = trees
burning = np.array([[0,0]])

# choose 5 randomly that are burning
rand = np.random.randint(N, size=5)
for i in range(len(rand)):
    tree = unburnt[i,:]
    print(tree)
    burning = np.concatenate((burning, np.array([tree])))
    unburnt = np.delete(unburnt, i, 0)

burning = np.delete(burning, 0, 0)
#%% Test plot forest 
plt.scatter(unburnt[:,0], unburnt[:,1], c = 'g', s = 1)
plt.scatter(burning[:,0], burning[:,1], c = 'r', s = 1)

# %%
# computes the distance with the 2-norm (squared) and initialization
def dist(v, u):
    sum = 0
    for i in range(len(v)):
        sum += (v[i]-u[i])**2
    return sum

#%% Test compute heat map and oxygen map
def update_heatmap_oxygenmap(unburnt, burning):
    temperature = np.zeros(omega)
    oxygen = np.zeros(omega)
    for i in range(temperature.shape[0]):
        for j in range(temperature.shape[1]):
            temp = 0
            oxy = 1
            for k in range(burning.shape[0]):
                temp += 600*np.exp(-1*dist(burning[k],[i,j])/8)
                oxy *= (1 - 0.5*np.exp(-1*dist(burning[k],[i,j]))/32)
            temperature[i,j] = temp
            oxygen[i,j] = oxy

    return [temperature, oxygen]

maps = update_heatmap_oxygenmap(unburnt,burning)
temperature = maps[0]
oxygen = maps[1]
#%% Test Heat Map
heat_map = plt.pcolormesh(temperature)
plt.colorbar(heat_map)
plt.title('Heat Map (degrees Celcius)')
#%% Test Heat Map 
oxygen_map = plt.pcolormesh(oxygen)
plt.colorbar(oxygen_map)
plt.title('Oxygen Map')
#%% Test adding trees 
def update_trees(unburnt, burning,temperature,oxygen):
    new_burning = burning
    new_unburnt = unburnt

    # update burning to burnt
    count = 0
    for i in range(len(burning)):
        if np.random.uniform(0,1) < 0.7:
            new_burning = np.delete(new_burning, i - count, 0)
            count += 1

    # update unburnt to burning
    count = 0
    for i in range(len(unburnt)):
        # get the coordinates of the tree 
        x = math.floor(unburnt[i,0])
        y = math.floor(unburnt[i,1])
        if temperature[x,y] > 300 and oxygen[x,y] > 0.2:
            tree = unburnt[i,:]
            new_burning = np.concatenate((new_burning, np.array([tree])))
            new_unburnt = np.delete(new_unburnt, i - count, 0)
            count += 1

    return [new_unburnt, new_burning]

new_trees = update_trees(unburnt,burning,temperature,oxygen)
unburnt = new_trees[0]
burning = new_trees[1]
#%% Test adding trees
plt.scatter(unburnt[:,0], unburnt[:,1], c = 'g', s = 1)
plt.scatter(burning[:,0], burning[:,1], c = 'r', s = 1)
# %%
time_step = 4
record_times = [1, 2, 3]

fig, ax = plt.subplots(len(record_times),3)

counter = 0
for t in range(time_step):
    print(t)
    new_maps = update_heatmap_oxygenmap(unburnt,burning)
    temperature = new_maps[0]
    oxygen = new_maps[1]
    new_trees = update_trees(unburnt,burning,temperature,oxygen)
    unburnt = new_trees[0]
    burning = new_trees[1]


    if t == record_times[counter]:
        ax[counter,0].scatter(unburnt[:,0], unburnt[:,1], c = 'g', s = 1)
        ax[counter,0].scatter(burning[:,0], burning[:,1], c = 'r', s = 1)
        #ax[counter,0].set_title('Forest Map at time t =' + str(record_times[counter]))


        ax[counter,1].pcolormesh(temperature)
        #ax[counter,1].set_title('Heat Map at time t = ' + str(record_times[counter]))

        ax[counter,2].pcolormesh(oxygen)
        #ax[counter,2].set_title('Oxygen Map at time t = ' + str(record_times[counter]))

        counter += 1

    

    
# %%
