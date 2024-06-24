# SCRIPT: Spatially varying mean velocity field
# MASTER STUDENT: Alberto Fontebasso
# PROJECT: Master Thesis
# SUPERVISOR: PD Dr. Meyer-Massetti
# LAB: IFD

import numpy as np
import matplotlib.pyplot as plt
import math
import random 
import sys
from scipy.interpolate import griddata
import multiprocessing

# NUMBER OF CORES
number_cores = 4
# READING FILE
points = np.loadtxt("Coordinates.txt", delimiter=',')
#st_dev_Y_values = np.loadtxt("st.txt", delimiter=',')
U1_values = np.loadtxt("U1.txt", delimiter=',')
U2_values = np.loadtxt("U2.txt", delimiter=',')
#value = griddata(points, st_dev_matrix, (-70, -60), method='cubic')

#------------------------------------------------------------------------------
# LOG CONDUCTIVITY PARAMETERS
var_Y = 1/16
st_dev_Y = np.sqrt(var_Y)
l_Y = 1
#------------------------------------------------------------------------------
# REAL LOOP FOR THE DETERMINATION OF THE TEMPORAL EVOLUTION VARIANCES
# Definition of initial and ending non-dimensional instants
t_prime_0 = 0 
t_prime_end = 40

# Time and time-step
length = 401
t_prime = np.linspace(t_prime_0,t_prime_end,length)
dt_prime = np.mean(np.diff(t_prime))

# Stochastic differential equation theta
# Definition of the drift coefficients
a_1_theta = 16/15
a_2_theta = (0.361**2)
# Definition of the diffusion coefficients
b_1_theta = 2/(np.sqrt(15))
# Stochastic coefficients theta
# Drift coefficient
drift_theta = lambda theta, y_prime: -(a_1_theta*theta + a_2_theta*y_prime)
# Diffusion coefficient
diffusion_theta = lambda st_dev_Y:  b_1_theta*st_dev_Y


# Stochastic differential equation v
# Definition of the drift coefficients
mu = lambda st_dev_Y: 0.34*(st_dev_Y**(2/1.9))
sigma = lambda st_dev_Y: 0.77*(st_dev_Y**(2/1.8))
alpha = lambda st_dev_Y: -1.2*(st_dev_Y**(2/5))
B = lambda st_dev_Y: 0.48*(st_dev_Y**(2/1.8)) #Diffusion coefficient
# Stochastic coefficients theta
# Drift coefficient
drift_v = lambda mu, sigma, alpha, B, v: (B**2)/(2*sigma)*((mu-v)/sigma + np.sqrt(2/math.pi)*alpha*math.exp(-(alpha**2)*((v-mu)**2)/(2*sigma**2))/math.erfc(alpha*(mu-v)/(np.sqrt(2)*sigma))) 
## Diffusion coefficient
diffusion_v = lambda B: B

# Solve Stochastic Differential Equation
n_trials = 5000

def stationary(k): 
    theta = 0
    v = 0
    y_prime = 0

    for i in range(1,length):  
        # Wiener increment 
        noise_theta = random.gauss(0.0, np.sqrt(dt_prime))
        noise_v = random.gauss(0.0, np.sqrt(dt_prime))
        
        theta += drift_theta(theta, y_prime)*dt_prime + diffusion_theta(st_dev_Y)*noise_theta    
        v += drift_v(mu(st_dev_Y), sigma(st_dev_Y), alpha(st_dev_Y), B(st_dev_Y), v)*dt_prime + diffusion_v( B(st_dev_Y) )*noise_v
        y_prime += np.exp(v)*(theta)*dt_prime
        
    return theta, v, y_prime

pool = multiprocessing.Pool(number_cores)
realizations_theta, realizations_v, realizations_y_prime = zip(*pool.map(stationary, range(n_trials)))

realizations_theta = np.asarray(realizations_theta)
realizations_v = np.asarray(realizations_v) 
realizations_y_prime = np.asarray(realizations_y_prime)        
#------------------------------------------------------------------------------     
# EXTRACTION OF THE STATIONARY PDF OF theta AND v
hist_theta, bins_theta = np.histogram(realizations_theta, bins=100)
bin_midpoints_theta = bins_theta[:-1] + np.diff(bins_theta)/2
cdf_theta = np.cumsum(hist_theta)
cdf_theta = cdf_theta / cdf_theta[-1]

hist_y_prime, bins_y_prime = np.histogram(realizations_y_prime, bins=100)
bin_midpoints_y_prime = bins_y_prime[:-1] + np.diff(bins_y_prime)/2
cdf_y_prime = np.cumsum(hist_y_prime)
cdf_y_prime = cdf_y_prime / cdf_y_prime[-1]

hist_v, bins_v = np.histogram(realizations_v, bins=100)
bin_midpoints_v = bins_v[:-1] + np.diff(bins_v)/2
cdf_v = np.cumsum(hist_v)
cdf_v = cdf_v / cdf_v[-1]
#------------------------------------------------------------------------------
hist_theta = []
bins_theta = []
hist_y_prime = []
bins_y_prime = []
hist_v = []
bins_v = []
realizations_theta = []
realizations_v = []
realizations_y_prime = []
# REAL LOOP FOR THE DETERMINATION OF THE TEMPORAL EVOLUTION VARIANCES
# Definition of initial and ending non-dimensional instants
t_prime_0 = 0 
t_prime_end = 100

# Time and time-step
length = 101
t_prime = np.linspace(t_prime_0,t_prime_end,length)
dt_prime = np.mean(np.diff(t_prime))

# Solve Stochastic Differential Equation
n_trials = 4000

def simulation(k):
    x1 = np.zeros(length) 
    x2 = np.zeros(length)
    velocity = np.zeros(length)

    # INITIAL POSITIONS
    x1[0] = q
    x2[0] = 15
        
    # INITIALIZATION FROM THE STATIONARY DISTRIBUTION
    prob_theta = np.random.rand()
    value_prob_theta = np.searchsorted(cdf_theta, prob_theta)
    theta = bin_midpoints_theta[value_prob_theta]
        
    prob_y_prime = np.random.rand()
    value_prob_y_prime = np.searchsorted(cdf_y_prime, prob_y_prime)
    y_prime = bin_midpoints_y_prime[value_prob_y_prime]
    
    prob_v = np.random.rand()
    value_prob_v = np.searchsorted(cdf_v, prob_v)
    v = bin_midpoints_v[value_prob_v]  
    velocity[0] = v
    
    for i in range(1,length):
        if (x1[i-1] > 98 or x2[i-1] > 98):
            break
        # EXTRACT STANDARD DEVIATION FROM POSITION
        U1 = griddata(points, U1_values, (x1[i-1], x2[i-1]), method='nearest')
        U2 = griddata(points, U2_values, (x1[i-1], x2[i-1]), method='nearest')
        U = np.sqrt(U1**2 + U2**2)
        
        if ( (U1 > 0 and U2 >= 0) or (U1 > 0 and U2 < 0) ):
            phi = np.arctan(U2/U1)
        elif ( (U1 < 0 and U2 >= 0) or (U1 < 0 and U2 < 0) ):
            phi = math.pi + np.arctan(U2/U1)
        elif (U1 == 0 and U2 >= 0):
            phi = math.pi/2
        elif (U1 == 0 and U2 < 0):
            phi = -math.pi/2
            
        # Wiener increment 
        noise_theta = random.gauss(0.0, np.sqrt(dt_prime))
        noise_v = random.gauss(0.0, np.sqrt(dt_prime))

        theta += drift_theta(theta, y_prime)*dt_prime + diffusion_theta(st_dev_Y)*noise_theta
        v += drift_v( mu(st_dev_Y), sigma(st_dev_Y) , alpha(st_dev_Y) , B(st_dev_Y), v)*dt_prime + diffusion_v(B(st_dev_Y) )*noise_v
        y_prime += np.exp(v)*theta*dt_prime
                
        x1[i] = x1[i-1] + np.exp(v)*np.cos(theta + phi)*dt_prime*l_Y
        x2[i] = x2[i-1] + np.exp(v)*np.sin(theta + phi)*dt_prime*l_Y       
        velocity[i] = np.exp(v)*U
        
    return x1, x2, velocity

vector = np.linspace(75,15,13)
count = 0
for q in vector:
    count += 1
    pool = multiprocessing.Pool(number_cores)
    realizations_x1, realizations_x2, realizations_v = zip(*pool.map(simulation, range(n_trials)))
        
    a = str(count)
    np.savetxt('{0}_x1.txt'.format(a), realizations_x1, delimiter=',')
    np.savetxt('{0}_x2.txt'.format(a), realizations_x2, delimiter=',')
    np.savetxt('{0}_velocity.txt'.format(a), realizations_v, delimiter=',')
#------------------------------------------------------------------------------