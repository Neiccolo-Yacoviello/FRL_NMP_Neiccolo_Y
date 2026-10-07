from integrator_constants import *
import numpy as np
from scipy.integrate import solve_ivp
from scipy.integrate import quad

"""Functions for the equations of motion in a 3DOF rocket model"""
# Returns the acceleration due to gravity at:
# height "h" (in meters) above the Earth's surface (for it's mean radius)
def gravity(h):
    g = g_Earth_S * (R_Earth / (R_Earth + h))**2
    return g

# Returns the thrust of the rocket engine at:
# time "t" (in seconds) where after burnout (t > 5 seconds) the thrust is zero
def thrust(t):
    if t < 5:
       T = thrust_LUT[t] # Retrieve thrust from LUT at time = t
    if t >= 5:
       return 0 # Burnout

# Returns the mass of the rocket at:
# time "t" (in seconds)
def mass(t):
    mass = m0 - quad(lambda t: -thrust(t) / (g_Earth_S * I_sp), 0, t)[0] # Integrate the mass flow rate from 0 to t to get the total change in mass at time t
    return mass

# Returns the magnitude of the drag force in the inertial frame at:
# velocity "v" (in m/s)
# height "h" (in meters)
# reference area "S" (in m^2)
def drag(v, h, S):
    drag = 0.5 * atmospheric_LUT[h][1] * v**2 * aerodynamic_LUT[v][1] * S
    return drag
    
# Returns the magnitude of the lift force in the inertial frame at:
# velocity "v" (in m/s)
# height "h" (in meters)
# reference area "S" (in m^2)
def lift(v, h, S):
    lift = 0.5 * atmospheric_LUT[h][1] * v**2 * aerodynamic_LUT[v][2] * S
    return lift

# Returns the weight of the rocket at:
# height "h" (in meters)
# time "t" (in seconds)
def weight(t, h):
    weight = mass(t) * gravity(h)
    return weight

# Returns the net force in the z-direction of the inertial frame at:
# time "t" (in seconds) 
# velocity "v" (in m/s)
# height "h" (in meters)
# reference area "S" (in m^2) 
# pitch angle "theta" (in radians measured from the vertical)
# angle of the velocity vector "upsilon" (in radians measured from the vertical)
def eom_F_z(t, h, v, S, theta, upsilon):
    F_z = thrust(t) * np.cos(theta) - weight(t, h) - drag(v, h, S) * np.cos(upsilon) - lift(v, h, S) * np.sin(upsilon)
    return F_z

# Returns the net force in the y-direction of the inertial frame at:
# time "t" (in seconds) 
# velocity "v" (in m/s)
# height "h" (in meters)
# reference area "S" (in m^2) 
# pitch angle "theta" (in radians measured from the vertical)
# angle of the velocity vector "upsilon" (in radians measured from the vertical)
def eom_F_y(t, h, v, S, theta, upsilon):
    F_y = thrust(t) * np.sin(theta) - drag(v, h, S) * np.sin(upsilon) + lift(v, h, S) * np.cos(upsilon)
    return F_y

# def eom_F_pitch():


# def model():
     

# def solve_equation():

# def get_derivations(eq, trial_point, step_size, num_steps):
#     klist = []
#     for i in range(1,6):
#         new_k = 