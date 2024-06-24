# SCRIPT: PMVP particle generator
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


points_X1 = np.loadtxt("Points_X1.txt", delimiter=',')
points_X2 = np.loadtxt("Points_X2.txt", delimiter=',')
U1_values = np.loadtxt("V_X1.txt", delimiter=',')
U2_values = np.loadtxt("V_X2.txt", delimiter=',')

vector = np.array([-3, -1, 0, 1, 3])

T_positions = np.array([0, 0])

def Autocovariance(x,y):
    st = 2.1366*np.exp( -np.sqrt( ( np.linalg.norm(x-y)/1 )**2 ) )
    return st

C_prime = Autocovariance(T_positions, T_positions)
C = np.zeros([2,2])
C[1][1] = C_prime


# NUMBER OF CORES
number_cores = 4

#LOG CONDUCTIVITY PARAMETERS
l_Y = 1

# DEFINITION OF THE DISCRETIZATION FOR THE PDF OF THE log-velocity
n_points = 50000;
v_values = np.linspace(-20,20,n_points)
dv = v_values[1] - v_values[0]
PDF_v = np.zeros(n_points)

# DEFINITION OF STARING AND ENDING TEMPORAL INSTANTS
t_prime_0 = 0 
t_prime_end = 30

# DISCRETIZED TIMELINE AND TIME-STEP SIZE
length = 120
t_prime = np.linspace(t_prime_0,t_prime_end,length)
dt_prime = np.mean(np.diff(t_prime))

# STOCHASTIC DIFFERENTIAL EQUATION FOR THE ANGLE PROCESS
# DEFINITION OF THE DRIFT COEFFICIENT CONSTANTS
a_1_theta = 16/15
a_2_theta = (0.361**2)
# DEFINITION OF THE DIFFUSION COEFFICIENT CONSTANTS
b_1_theta = 2/(np.sqrt(15))
# DRIFT COEFFICIENT
drift_theta = lambda theta, y_prime: -(a_1_theta*theta + a_2_theta*y_prime)
# DIFFUSION COEFFICIENT
diffusion_theta = lambda st_dev_Y:  b_1_theta*st_dev_Y

# STOCHASTIC DIFFERENTIAL EQUATION FOR THE LOG-VELOCITY PROCESS
# DEFINITION OF THE DRIFT COEFFICIENTS PARAMETERS
mu = lambda st_dev_Y: 0.34*(st_dev_Y**(2/1.9))
sigma = lambda st_dev_Y: 0.77*(st_dev_Y**(2/1.8))
alpha = lambda st_dev_Y: -1.2*(st_dev_Y**(2/5))
# DEFINITION OF THE DIFFUSION COEFFICIENTS PARAMENTERS
B = lambda st_dev_Y: 0.48*(st_dev_Y**(2/1.8))
# DRIFT COEFFICIENT
drift_v = lambda mu, sigma, alpha, B, v: (B**2)/(2*sigma)*((mu-v)/sigma + np.sqrt(2/math.pi)*alpha*math.exp(-(alpha**2)*((v-mu)**2)/(2*sigma**2))/math.erfc(alpha*(mu-v)/(np.sqrt(2)*sigma))) 
# DIFFUSION COEFFICIENT
diffusion_v = lambda B: B

# DEFINITION OF THE NUMBER OF REALIZATIONS
n_trials = 100000

#------------------------------------------------------------------------------
def simulation(k):

    # PARTICLE PATH INITIALIZATION
    x1 = np.zeros(length) 
    x2 = np.zeros(length)
    velocity_x1 = np.zeros(length)
    velocity_x2 = np.zeros(length)
    
    # INITIAL POSITIONS
    x1[0] = -6
    x2[0] = vector[q]
    
    # INITIALIZATION FROM THE STATIONARY DISTRIBUTION
    theta = random.gauss(0.0, st_dev_theta)  
    y_prime = random.gauss(0.0, st_dev_y_prime)
    
    prob_v = np.random.rand()
    value_prob_v = np.searchsorted(cdf_v, prob_v)
    v = v_values[value_prob_v]
    velocity_x1[0] = np.exp(v)*np.cos(theta)
    velocity_x2[0] = np.exp(v)*np.sin(theta)
    
    for i in range(1,length):
        if ( x1[i-1] > 6 or x2[i-1] < -6 or x2[i-1] > 6):
            break
        
        C[0][1] = Autocovariance(np.array([x2[i-1], x1[i-1]]), T_positions)
        C[0][0] = Autocovariance(np.array([x2[i-1], x1[i-1]]), np.array([x2[i-1], x1[i-1]]))
        C[1][0] = C[0][1]
        st_dev_Y = np.sqrt( np.linalg.det(C)/C_prime )        
        
        
        U1 = griddata(points_X1, U1_values, (x2[i-1], x1[i-1]), method='nearest')
        U2 = griddata(points_X2, U2_values, (x2[i-1], x1[i-1]), method='nearest')
#        U1 = 1
#        U2 = 0
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
                
        x1[i] = x1[i-1] + np.exp(v)*np.cos(theta + phi)*dt_prime
        x2[i] = x2[i-1] + np.exp(v)*np.sin(theta + phi)*dt_prime
        velocity_x1[i] = np.exp(v)*np.cos(theta + phi)*U
        velocity_x2[i] = np.exp(v)*np.sin(theta + phi)*U
    
        
    return x1, x2, velocity_x1, velocity_x2
#------------------------------------------------------------------------------

for q in range(5):
    # EXTRACT STANDARD DEVIATION OF THE CONDUCTIVITY AT THE INJECTION POINTS
#    st_dev_Y = griddata(points, st_dev_Y_values, (injection_x1[q], injection_x2[q]), method='nearest')
    
    C[0][1] = Autocovariance(np.array([vector[q], -6]), T_positions)
    C[0][0] = Autocovariance(np.array([vector[q], -6]), np.array([vector[q], -6]))
    C[1][0] = C[0][1]
    st_dev_Y = np.sqrt( np.linalg.det(C)/C_prime )
    
    # DEFINE STANDARD DEVIATION FOR THE THETA AND Y_PRIME STATIONARY DISTRIBUTIONS
    st_dev_theta = 0.37*(st_dev_Y**(2/1.92))
    if (st_dev_Y**2 >= 0 and st_dev_Y**2 < 1 ):
        st_dev_y_prime = 1.0829*(st_dev_Y**(2*0.5442))
    elif (st_dev_Y**2 >= 1):
        st_dev_y_prime = 0.5943*st_dev_Y**2 + 0.5009
        
    for p in range(0,n_points):
        PDF_v[p] = 1/(np.sqrt(2*math.pi)*sigma(st_dev_Y))*np.exp(-((v_values[p] - mu(st_dev_Y))**2 )/(2*sigma(st_dev_Y)**2))*math.erfc(-alpha(st_dev_Y)*(v_values[p] - mu(st_dev_Y))/(np.sqrt(2)*sigma(st_dev_Y))) 
    cdf_v = np.cumsum(PDF_v*dv)
    cdf_v = cdf_v / cdf_v[-1]
        
    pool = multiprocessing.Pool(number_cores)
    realizations_x1, realizations_x2, realizations_velocity_x1, realizations_velocity_x2 = zip(*pool.map(simulation, range(n_trials)))
    
    realizations_x1 = np.asarray(realizations_x1)
    realizations_x2 = np.asarray(realizations_x2) 
    realizations_velocity_x1 = np.asarray(realizations_velocity_x1)
    realizations_velocity_x2 = np.asarray(realizations_velocity_x2)

    a = str(q)
    np.savetxt('./FILE/{0}_x1.txt'.format(a), realizations_x1, delimiter=',')
    np.savetxt('./FILE/{0}_x2.txt'.format(a), realizations_x2, delimiter=',')
    np.savetxt('./FILE/{0}_velocity_x1.txt'.format(a), realizations_velocity_x1, delimiter=',')
    np.savetxt('./FILE/{0}_velocity_x2.txt'.format(a), realizations_velocity_x2, delimiter=',')
