#%%
import numpy as np
import matplotlib.pyplot as plt


#%% 
# data points drawn from uniform(0,1)
v = np.zeros(500)
for i in range(len(v)):
    v[i] = np.random.uniform()

plt.hist(v, bins=50)
#%%
# this function takes in a value from [0,1] and output the 
# corresponding input from CDF of a Laplace distribution
def inverseLaplace(xi, mu, w):
    if xi < 0.5: 
        return mu + w*np.log(2*xi)
    # Laplace distribution is symmetric, so the CDF is symmetric
    # about 0.5 
    else: 
        return  mu - w*np.log(2*(1 - xi))

w = np.zeros(len(v))
for i in range(len(v)):
    w[i] = inverseLaplace(v[i], 1, 2)



# the model of the ER

# %% 
# Time-independent parameters
# The maximum rate for which the waiting room help people
V1 = 12

# How quickly the waiting room saturates
K1 = 20

# The maximum rate for which the lanes help people
V2 = 6

# How quickly the lanes saturate
K2 = 10

# distribution of resources between green and red lane
alpha = 0.5

# The maximum rate for which the hospital let people out
V3 = 5

# How quickly the hospital saturate
K3 = 120

# %%
# Special Functions 
# Box-car function 
def boxcar(t, start, end, amplitutde):
    if t <= end and t >= start:
        return amplitutde
    else: return 0

# Sine function truncated at start time, end time, and values below 0
def truncatedSine(t, start, end, middle, amplitutde):
    if amplitutde * np.sin(t) + middle >= 0 and t >= start and t <= end:
        return middle + amplitutde * np.sin(t)
    else: return 0

# %%
# Time-dependent parameters
# Flux from home to waiting room with non-critical condition
def flux_nc(t):
    return 2 + boxcar(t, 10, 140, 1)

# Flux from home to waiting room with critical condition
def flux_c(t):
    return 0.5 + boxcar(t, 10, 120, 1)

# Flux from waiting room to red lane 
def flux_wr(W1, W2, t):
    W = W1 + W2
    return V1*(W1/W)*(W/(W + K1))

# Flux from waiting room to green lane
def flux_wg(W1, W2, t):
    W = W1 + W2
    return V1*(W2/W)*(W/(W + K1))

# Flux to red lane from ambulances
def flux_am(t):
    return 0.2 + truncatedSine(t, 0, 120, 0.2, 4)

# Flux from red lane to hospital
def flux_rh(R, G, t):
    P = R + G
    return V2*((R/P)**alpha)*(P/(P+K2)) 

# Flux from green lane back home
def flux_gh(R, G, t):
    P = R + G
    return V2*(1 - (R/P)**alpha)*(P/(P+K2)) 

# Flux from hospital back home
def flux_hh(H, t):
    return V3*(H/(H+K3))
#%%
# Model
from scipy.integrate import odeint
import numpy as np

# C[0] = Waiting room (critical), C[1] = Waiting room (non-critical),
# C[1] = Red Lane, C[2] = Green Lane, C[3] = Hospital
def model(C, t):
    waiting_room_c = flux_c(t) - flux_wr(C[0],C[1],t)
    waiting_room_nc = flux_nc(t) - flux_wg(C[0],C[1],t)
    red_lane = flux_wr(C[0],C[1],t) + flux_am(t) - flux_rh(C[2], C[3], t)
    green_lane = flux_wg(C[0],C[1],t) - flux_gh(C[2], C[3], t)
    hospital = flux_rh(C[2], C[3], t) - flux_hh(C[4], t)
    return [waiting_room_c, waiting_room_nc, red_lane, green_lane, hospital]

time = np.linspace(0, 216, 1001)

initial_cond = [1, 0, 1, 1, 40]

output = odeint(model, initial_cond, time)
fig, ax = plt.subplots()
hs = ax.twinx()

ax.plot(time, output[:,0] + output[:,1], label = 'Waiting Room')
ax.plot(time, output[:,2], 'r', label = 'Red Lane')
ax.plot(time, output[:,3], 'g', label = 'Green Lane')
hs.plot(time, output[:,4], 'y', label = 'Hospital')
ax.legend(loc = 'upper left')
hs.legend(loc = 'upper right')
ax.set_ylabel('number of patients in ER')
ax.set_ylim(0, 100)
hs.set_ylabel('number of patients in hospital')
hs.set_ylim(0, 100)
fig.suptitle('ER and hospital during a storm')

# %%
# Stochastic Model

nstep = 2000
stime = np.zeros(nstep)
stime[0] = 0

# stores: 0. Waiting room (critical), 1. Waiting room(non-critical),
# 2. Red Lane, 3. Green Lane, 4. Hospital
X = np.zeros([5,nstep])
X[:,0] = initial_cond

# Stoichiometric Matrix
S = np.array([[-1,0,0,0,0],[0,-1,0,0,0],[1,0,-1,0,0],[0,1,0,-1,0],[0,0,1,0,-1]])

# propensity as the out-flux of each compartment
def propensity(X, t):
    waiting_room_c = flux_wr(X[0],X[1],t)
    waiting_room_nc = flux_wg(X[0],X[1],t)
    red_lane = flux_rh(X[2],X[3],t)
    green_lane = flux_gh(X[2],X[3],t)
    hospital = flux_hh(X[4],t)
    return [waiting_room_c, waiting_room_nc, red_lane, green_lane, hospital]

# %%
def checkProp(P, X):
    for j in range(len(P)):
        if np.isnan(P[j]) or P[j] < 0 or X[j] <= 0:
            P[j] = 0


for i in range(nstep-1):
    # Draw the next reaction time according to total propensity
    prop = propensity(X[:,i],stime[i]) 
    checkProp(prop, X[:,i])
    prop0 = sum(prop)
    xi = np.random.uniform(0,1)
    incr = -1/prop0 * np.log(1 - xi)
    stime[i+1] = stime[i] + incr

    # Draw next reaction proportional to each propensity
    xi = np.random.uniform(0,1)
    count = 0 
    xi -= prop[0]/prop0
    while (xi > 0):
        count += 1
        xi -= prop[count]/prop0
    #print(count)
    
    # Update the number of people
    X[:,i+1] = X[:,i] + S[:,count]

    # New arrivals (assume homogeneous Poisson Process)
    X[0,i+1] += np.random.poisson(flux_c(stime[i])*incr)
    X[1,i+1] += np.random.poisson(flux_nc(stime[i])*incr)
    X[2,i+1] += np.random.poisson(flux_am(stime[i])*incr)

# %% 
# plots the outcome of the stochastic model
plt.plot(stime,X[0,:] + X[1,:], label = 'Waiting Room')
plt.plot(stime,X[2,:], 'r', label = 'Red Lane')
plt.plot(stime,X[3,:], 'g', label = 'Green Lane')
plt.ylim(0,120)
plt.plot(stime,X[4,:], 'y', label = 'Hospital')
plt.title('ER and Hospital with storm for the first 5 days')
plt.legend()
