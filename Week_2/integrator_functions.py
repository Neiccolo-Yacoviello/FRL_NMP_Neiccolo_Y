from integrator_constants import *
import numpy as np
from scipy.integrate import quad

"""Functions for the inputs of the equations of motion in a 3DOF rocket model"""
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
       T = 0 # Burnout (later will be entirely controlled by LUT)
    return T

# Returns the mass flow rate of the rocket at:
# time "t" (in seconds)
def mass_flow_rate(t):
    m_dot = -thrust(t) / (g_Earth_S * I_sp)
    return m_dot

# Returns the mass of the rocket at:
# time "t" (in seconds)
def mass(t):
    mass = m0 + quad(lambda x: mass_flow_rate(x), 0, t)[0] # Integrate the mass flow rate from 0 to t to get the total change in mass at time t
    return mass

# Returns the magnitude of the drag force in the inertial frame at:
# height "h" (in meters)
# velocity "v" (in m/s)
# reference area "S" (in m^2)
def drag(h, v):
    drag = 0.5 * atmospheric_LUT[h][1] * v**2 * aerodynamic_LUT[v][1] * S_area
    return drag
    
# Returns the magnitude of the lift force in the inertial frame at:
# height "h" (in meters)
# velocity "v" (in m/s)
# reference area "S" (in m^2)
def lift(h, v):
    lift = 0.5 * atmospheric_LUT[h][1] * v**2 * aerodynamic_LUT[v][2] * S_area
    return lift

# Returns the weight of the rocket at:
# height "h" (in meters)
# time "t" (in seconds)
def weight(t, h):
    weight = mass(t) * gravity(h)
    return weight

"""EOMs (Equations of Motion) that describe each degree of freedom in a 3DOF rocket model"""
# Returns the acceleration in the z-direction of the inertial frame at:
# time "t" (in seconds)
# height "h" (in meters)
# velocity "v" (in m/s)
# reference area "S" (in m^2) 
# pitch angle "theta" (in radians measured from the vertical)
# angle of the velocity vector "upsilon" (in radians measured from the vertical)
def eom_a_z(t, h, v, theta, upsilon):
    a_z = (thrust(t) * np.cos(theta) - weight(t, h) - drag(h, v) * np.cos(upsilon) - lift(h, v) * np.sin(upsilon))/mass(t)
    return a_z

# Returns the acceleration in the y-direction of the inertial frame at:
# time "t" (in seconds)
# height "h" (in meters)
# velocity "v" (in m/s)
# reference area "S" (in m^2) 
# pitch angle "theta" (in radians measured from the vertical)
# angle of the velocity vector "upsilon" (in radians measured from the vertical)
def eom_a_y(t, h, v, theta, upsilon):
    a_y = (thrust(t) * np.sin(theta) - drag(h, v) * np.sin(upsilon) + lift(h, v) * np.cos(upsilon))/mass(t)
    return a_y

# Returns the acceleration of the pitch angle of the inertial frame at:
# height "h" (in meters)
# velocity "v" (in m/s)
# reference area "S" (in m^2)
# change in pitch angle "delta_theta" (in radians)
def eom_a_pitch(h, v, delta_theta):
    a_pitch = (aerodynamic_LUT[v][2] * (aerodynamic_LUT[v][3] - aerodynamic_LUT[v][1]) + 0.25 * atmospheric_LUT[h][1] * S_area * d**2 * C_m_q * delta_theta)/inertia_tensor[0][0]
    return a_pitch

"""Functions that describe each RK45 step"""
# Returns the derivative of the state vector at the given time and state for the RK45 integrator
# time "t" (in seconds)
# state vector "x" containing [z, y, v_z, v_y, theta, theta_dot]
def f(t, x):
    y, z, v_y, v_z, theta, theta_dot = x

    v = np.sqrt(v_z**2 + v_y**2)

    a_z = eom_a_z(t, z, v, theta, np.arctan2(v_y, v_z))
    a_y = eom_a_y(t, z, v, theta, np.arctan2(v_y, v_z))
    a_pitch = eom_a_pitch(z, v, theta)

    return np.array([v_y, v_z, a_y, a_z, theta_dot, a_pitch, mass_flow_rate(t)])

# Returns the list of derivative approximations at the given trial point for:
# the equation of motion "eq" (either a_z, a_y, or a_pitch)
# the time at which the derivative is evaluated "trial_point" (in seconds)
# the step size for the finite difference approximation "step_size" (in seconds)
def rk_45_step(trial_point, step_size, state_vector):
    k1 = step_size * f(trial_point + step_size * b_tableau[0][0], state_vector)
    k2 = step_size * f(trial_point + step_size * b_tableau[1][0], state_vector + k1 * b_tableau[1][1])
    k3 = step_size * f(trial_point + step_size * b_tableau[2][0], state_vector + k1 * b_tableau[2][1] + k2 * b_tableau[2][2])
    k4 = step_size * f(trial_point + step_size * b_tableau[3][0], state_vector + k1 * b_tableau[3][1] + k2 * b_tableau[3][2] + k3 * b_tableau[3][3])
    k5 = step_size * f(trial_point + step_size * b_tableau[4][0], state_vector + k1 * b_tableau[4][1] + k2 * b_tableau[4][2] + k3 * b_tableau[4][3] + k4 * b_tableau[4][4])
    k6 = step_size * f(trial_point + step_size * b_tableau[5][0], state_vector + k1 * b_tableau[5][1] + k2 * b_tableau[5][2] + k3 * b_tableau[5][3] + k4 * b_tableau[5][4] + k5 * b_tableau[5][5])

    x4 = state_vector + k1 * b_tableau[6][1] + k2 * b_tableau[6][2] + k3 * b_tableau[6][3] + k4 * b_tableau[6][4] + k5 * b_tableau[6][5]
    x5 = state_vector + k1 * b_tableau[7][1] + k2 * b_tableau[7][2] + k3 * b_tableau[7][3] + k4 * b_tableau[7][4] + k5 * b_tableau[7][5] + k6 * b_tableau[7][6]
    return x5, x5 - x4

# def check_step():
