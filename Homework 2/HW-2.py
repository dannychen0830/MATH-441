import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Problem 1 - Rise and Fall of Dynasties

# vector of parameters
#    a, b,    c,   d,    e,   f,   g,     h,    K, m,   r
#    0, 1,    2,   3,    4,   5,   6,     7,    8, 9,   10
k = [1, 0.17, 0.4, 0.42, 1.2, 0.1, 0.009, 0.1 , 1, 0.4, 1]

# model:F = p[0], B = p[1], R = p[2]
# First model for problem 1, where the equations are from the ones given
def P1_1(p, t, k):
    dFdt = k[10]*p[0]*(1 - p[0]/k[8]) - k[1]*p[0]*p[1]/(p[0] + k[2]) - k[7]*p[0]*p[2]
    dBdt = k[4]*k[1]**p[0]*p[1]/(p[0] + k[2]) - k[9]*p[1] - k[2]*p[1]*p[2]/(k[3] + p[2])
    dRdt = k[5]*k[1]**p[0]*p[1]/(p[0] + k[2]) - k[6]*p[2]
    return [dFdt, dBdt, dRdt]

# Second model for problem 1, now the rulers benefit from taxing
def P1_2(p, t, k):
    dFdt = k[10]*p[0]*(1 - p[0]/k[8]) - k[1]*p[0]*p[1]/(p[0] + k[2]) - k[7]*p[0]*p[2]
    dBdt = k[4]*k[1]**p[0]*p[1]/(p[0] + k[2]) - k[9]*p[1] - k[2]*p[1]*p[2]/(k[3] + p[2])
    dRdt = k[5]*k[1]**p[0]*p[1]/(p[0] + k[2]) - k[6]*p[2] + k[7]*p[0]*p[2]
    return [dFdt, dBdt, dRdt]

# initial values
P1_0 = [0.5, 0.2, 0.01] 

# time increments
t_1 = np.linspace(0, 10, 101)

# Plot solution function
def plotSol_1(obj, func):
    sol1 = odeint(func,P1_0, t_1, args=(k,))
    obj.plot(t_1, sol1[:,0],'b', label = 'Farmers')
    obj.plot(t_1, sol1[:,1],'r', label = 'Bandits')
    obj.plot(t_1, sol1[:,2],'g', label = 'Rulers')
    obj.legend(loc = 'upper left')
    obj.set_xlabel("time")
    obj.set_ylabel("normalized population ")   
    obj.set_ylim(0,1.05)

# Testing the effect of different parameters:
# Test number 1: Larger growth in Bandits
fig3, axs3 = plt.subplots(3, 3, constrained_layout = True)
fig3.suptitle('Effect of Parameters on Farmer-Bandit Dynamics')
banditRate = k[4]
farmerRate = k[0]
k[4] -= 1
k[0] -= 0.5
for i in range(3):
    k[0] = farmerRate - 0.5
    for j in range(3):
        plotSol_1(axs3[i,j],P1_1)
        axs3[i,j].set_title("a = " + str(round(k[0],2)) + ", e = " + str(round(k[4],2)))
        axs3[i,j].get_legend().remove()
        axs3[i,j].set_ylabel("population")
        k[0] += 0.5
    k[4] += 1
k[4] = banditRate
k[0] = farmerRate


# Comparision between not / gaining advantage from tax
fig1, axs1 = plt.subplots(1,2, constrained_layout = True)

plotSol_1(axs1[0], P1_1)
axs1[0].set_title("Population Dynamics when Rulers \n don't Benefit from Tax")

plotSol_1(axs1[1], P1_2)
axs1[1].set_title("Population Dynamics when Rulers \n Benefit from Tax")


# Test the effect of amount of taxes on the dynamics of the system

fig2, axs2 = plt.subplots(nrows= 2,ncols= 3,constrained_layout = True)
fig2.suptitle('Population Dynamics with Varying Amount of Taxation')
for i in range(2):
    for j in range(3):
        if i == 1 and j == 2:
            axs2[i,j].axis('off')
            break
        sol = odeint(P1_2, P1_0, t_1, args=(k,))
        axs2[i,j].plot(t_1, sol[:,0],'b', label = 'Farmers')
        axs2[i,j].plot(t_1, sol[:,1],'r', label = 'Bandits')
        axs2[i,j].plot(t_1, sol[:,2],'g', label = 'Rulers')
        axs2[i,j].set_xlabel('time')
        axs2[i,j].set_ylabel('normalized population')
        axs2[i,j].set_title('h = ' + str(round(k[7],2)))
        axs2[i,j].set_ylim(0,1.05)
        k[7] += 0.1
lines, labels = fig2.axes[0].get_legend_handles_labels()
fig2.legend(lines, labels, loc = 'lower right',bbox_to_anchor=(0.9,0.15))
k[7] = 0.1

# Test the effect of the rulers effectively chasing down bandits
fig4, axs4 = plt.subplots(1,3, constrained_layout = True)
fig4.suptitle("The Effect of Ruler's Attitude on Bandits")
k[2] -= 0.3
for i in range(3):
    plotSol_1(axs4[i], P1_2)
    axs4[i].set_title("c = " + str(round(k[2],2)))
    k[2] += 0.3


# -----------
# Problem 4 - Age Structured Population Model

# Mortality (mu)
def mu(j):
    return np.exp(j/16.07) / (199 + np.exp(j/16.07))

# Feunidity (beta)
def beta(j):
    return ((j/0.75)**40) * np.exp(-j/0.75) / 3.871318360247655e+47

def EstimateB():
    sum = 0
    prod = 1
    for j in range(1, 121):
        for k in range(j):
            prod *= 1 - mu(k)
        sum += beta(j)* prod
        print(j,':', beta(j))
        prod = 1

def LeslieMatrix():
    L = np.zeros((120,120))
    for i in range(120):
        L[0,i] = beta(i+1)
    for j in range(1,120):
        L[j, j-1] = 1 - mu(j) 
    return L

# Create an random population distribution and normalize it
P4_0 = np.random.uniform(0,1,120)
P4_0 = P4_0 / np.sum(P4_0)

# The age axis 
age = np.linspace(0,119,120)

# times we are testing
time = [0, 1, 20, 40, 60, 120, 240, 500]

L = LeslieMatrix()

fig3, axs3 = plt.subplots(2,4, constrained_layout = True)
fig3.suptitle("Human Population Compartmentalized by Age over time")
for i in range(8):
    R = np.linalg.matrix_power(L,time[i])
    pop = np.dot(R,P4_0)
    axs3[int(i/4), int(i%4)].plot(age, pop / np.sum(pop))
    axs3[int(i/4), int(i%4)].set_xlabel("age")
    axs3[int(i/4), int(i%4)].set_ylabel("noralized population")
    axs3[int(i/4), int(i%4)].set_title("After " + str(time[i]) + " years")
    axs3[int(i/4), int(i%4)].set_ylim([0, 0.05])
    axs3[int(i/4), int(i%4)].set_xlim([0,120])

plt.show()




