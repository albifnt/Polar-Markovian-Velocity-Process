# SCRIPT: v-theta RW model implied by perturbation theory
# MASTER STUDENT: Alberto Fontebasso
# PROJECT: Master Thesis
# SUPERVISOR: PD Dr. Meyer-Massetti
# LAB: IFD

import numpy as np
import matplotlib.pyplot as plt
import random 


# LOG CONDUCTIVITY PARAMETERS
var_Y = 4
st_dev = np.sqrt(var_Y)

# We do everything with adimensional variables
# Definition of initial and ending adimensional instants
t_prime_0 = 0 
t_prime_end = 40

# Time and time-step
length = 401
t_prime = np.linspace(t_prime_0,t_prime_end,length)
dt_prime = np.mean(np.diff(t_prime))

# Stochastic differential equation theta
# Definition of the drift coefficients
a_1_theta = 16/15
a_2_theta = 0.361**2

# Definition of the diffusion coefficients
b_1_theta = 2/(np.sqrt(15))*st_dev

# Stochastic differential equation v
# Definition of the drift coefficients
a_1_v = - 8/15
# Definition of the diffusion coefficients
b_1_v = np.sqrt(2/5)*st_dev


# Stochastic coefficients theta
# Drift coefficient
drift_theta = lambda theta,y_prime: - (a_1_theta*theta + a_2_theta*y_prime)
# Diffusion coefficient
diffusion_theta = b_1_theta

# Stochastic coefficients theta
# Drift coefficient
drift_v = lambda v: a_1_v*v
# Diffusion coefficient
diffusion_v = b_1_v


# Solve Stochastic Differential Equation
n_trials = 5000
realizations_theta = np.zeros((n_trials,1))
realizations_v = np.zeros((n_trials,1))
realizations_y_prime = np.zeros((n_trials,1))

for k in range(0,n_trials):
    theta = 0
    v = 0
    y_prime = 0

    for i in range(1,length):
        # Wiener increment 
        noise_theta = random.gauss(0.0, np.sqrt(dt_prime))
        noise_v = random.gauss(0.0, np.sqrt(dt_prime))
        
        theta += drift_theta(theta, y_prime)*dt_prime + diffusion_theta*noise_theta
        v += drift_v(v)*dt_prime + diffusion_v*noise_v
        y_prime += np.exp(v)*(theta)*dt_prime
        
        if (i == length-1):
            realizations_theta[k][0] = theta
            realizations_v[k][0] = v
            realizations_y_prime[k][0] = y_prime

#------------------------------------------------------------------------------     
# EXTRACTION OF THE STATIONARY PDF OF theta AND v
hist_theta, bins_theta = np.histogram(realizations_theta[:,0], bins=100)
bin_midpoints_theta = bins_theta[:-1] + np.diff(bins_theta)/2
cdf_theta = np.cumsum(hist_theta)
cdf_theta = cdf_theta / cdf_theta[-1]

hist_y_prime, bins_y_prime = np.histogram(realizations_y_prime[:,0], bins=100)
bin_midpoints_y_prime = bins_y_prime[:-1] + np.diff(bins_y_prime)/2
cdf_y_prime = np.cumsum(hist_y_prime)
cdf_y_prime = cdf_y_prime / cdf_y_prime[-1]

hist_v, bins_v = np.histogram(realizations_v[:,0], bins=100)
bin_midpoints_v = bins_v[:-1] + np.diff(bins_v)/2
cdf_v = np.cumsum(hist_v)
cdf_v = cdf_v / cdf_v[-1]
#------------------------------------------------------------------------------

# Definition of initial and ending adimensional instants
t_prime_0 = 0 
t_prime_end = 32

# Time and time-step
length = 513
t_prime = np.linspace(t_prime_0,t_prime_end,length)
dt_prime = np.mean(np.diff(t_prime))

# Solve Stochastic Differential Equation
n_trials = 5000
realizations_x1 = np.zeros((n_trials,length))
realizations_x2 = np.zeros((n_trials,length))

for k in range(0,n_trials):
    x1 = np.zeros(length) 
    x2 = np.zeros(length)
    
    prob_theta = np.random.rand()
    value_prob_theta = np.searchsorted(cdf_theta, prob_theta)
    theta = bin_midpoints_theta[value_prob_theta]
    
    prob_y_prime = np.random.rand()
    value_prob_y_prime = np.searchsorted(cdf_y_prime, prob_y_prime)
    y_prime = bin_midpoints_y_prime[value_prob_y_prime]
    
    prob_v = np.random.rand()
    value_prob_v = np.searchsorted(cdf_v, prob_v)
    v = bin_midpoints_v[value_prob_v]  

    for i in range(1,length):
        # Wiener increment 
        noise_theta = random.gauss(0.0, np.sqrt(dt_prime))
        noise_v = random.gauss(0.0, np.sqrt(dt_prime))
        
        theta += drift_theta(theta, y_prime)*dt_prime + diffusion_theta*noise_theta
        v += drift_v(v)*dt_prime + diffusion_v*noise_v
        y_prime += np.exp(v)*theta*dt_prime
        
        x1[i] = x1[i-1] + np.exp(v)*np.cos(theta)*dt_prime 
        x2[i] = x2[i-1] + np.exp(v)*np.sin(theta)*dt_prime
        
        realizations_x1[k][i] = x1[i]
        realizations_x2[k][i] = x2[i]

variance = np.zeros(length)
for z in range(0,length):
    variance[z] = np.var(realizations_x1[:,z], ddof=1)

#plt.plot(t_prime,variance)
np.savetxt('Variance_x_116.txt', variance, delimiter=',')
#variance = np.zeros(length)
#for z in range(0,length):
#    variance[z] = np.var(realizations_x2[:,z], ddof=1)

#np.savetxt('Variance_y_116.txt', variance, delimiter=',')
#plt.plot(t_prime,variance)
#np.savetxt('Time.txt', t_prime, delimiter=',')

