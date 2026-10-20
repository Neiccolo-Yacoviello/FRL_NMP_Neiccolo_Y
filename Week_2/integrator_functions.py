from integrator_constants import *
import numpy as np

"""Functions for the equations of motion in a 3DOF rocket model"""
# Returns the acceleration due to gravity at:
# height "h" (in meters) above the Earth's surface (for it's mean radius)
def gravity(h):
    g = g_Earth_S * (R_Earth / (R_Earth + h))**2
    return g

# Returns the thrust of the rocket engine at:
# time "t" (in seconds) where after burnout (t > 5 seconds) the thrust is zero
def thrust(t):
    T = thrust_LUT(t)
    return T

# Returns the mass flow rate of the rocket at:
# time "t" (in seconds)
def mass_flow_rate(t):
    m_dot = -thrust(t) / (g_Earth_S * I_sp)
    return m_dot

# Returns the dynamic pressure on the rocket at:
# height "h" (in meters)
# velocity "v" (in m/s)
def q(h, v):
    Q = 0.5 * atmospheric_LUT(h)[0] * v**2
    return Q

# Returns the mach number of the rocket at:
# height "h" (in meters)
# velocity "v" (in m/s)
def get_mach(h, v):
    temp = atmospheric_LUT(h)[2] + 273.15
    local_mach = np.sqrt(temp * air_specific_gas_constant * air_specific_heat_ratio)
    return np.abs(v)/local_mach

# Returns the coeffcient of lift (based on the Normal force coefficient and Axial force coefficient) at:
# time "t" (in seconds)
# height "h" (in meters)
# velocity "v" (in m/s)
# the pitch of the rocket "pitch"
# the angle of the velocity vector "upsilon"
def get_C_L(t, h, v, pitch, upsilon):
    AoA = pitch - upsilon
    mach = get_mach(h, v)
    C_L = aerodynamic_LUT(t, mach)[3] * np.cos(AoA) - aerodynamic_LUT(t, mach)[2] * np.sin(AoA)
    return C_L

# Returns the magnitude of the drag force in the inertial frame at:
# height "h" (in meters)
# velocity "v" (in m/s)
def drag(t, h, v):
    drag = q(h, v) * aerodynamic_LUT(t, get_mach(h, v))[1] * S_area
    return drag
    
# Returns the magnitude of the lift force in the inertial frame at:
# time "t" (in seconds)
# height "h" (in meters)
# velocity "v" (in m/s)
# the pitch of the rocket "pitch"
# the angle of the velocity vector "upsilon"
def lift(t, h, v, pitch, upsilon):
    lift = q(h, v) * get_C_L(t, h, v, pitch, upsilon) * S_area
    return float(lift)

# Returns the weight of the rocket at:
# height "h" (in meters)
# mass "m" (in kg)
def weight(h, m):
    weight = m * gravity(h)
    return weight

"""EOMs (Equations of Motion) that describe each degree of freedom in a 3DOF rocket model"""
# Returns the acceleration in the z-direction of the inertial frame at:
# time "t" (in seconds)
# height "h" (in meters)
# velocity "v" (in m/s)
# pitch angle "pitch" (in radians measured from the vertical)
# angle of the velocity vector "upsilon" (in radians measured from the vertical)
# current mass of the rocket "m" (in kg)
def eom_a_z(t, h, v, pitch, upsilon, m):
    a_z = (thrust(t) * np.cos(pitch) - weight(h, m) - drag(t, h, v) * np.cos(upsilon) - lift(t, h, v, pitch, upsilon) * np.sin(upsilon))/m
    return a_z

# Returns the acceleration in the y-direction of the inertial frame at:
# time "t" (in seconds)
# height "h" (in meters)
# velocity "v" (in m/s)
# pitch angle "pitch" (in radians measured from the vertical)
# angle of the velocity vector "upsilon" (in radians measured from the vertical)
# current mass of the rocket "m" (in kg)
def eom_a_y(t, h, v, pitch, upsilon, m):
    a_y = (thrust(t) * np.sin(pitch) - drag(t, h, v) * np.sin(upsilon) + lift(t, h, v, pitch, upsilon) * np.cos(upsilon))/m
    return a_y

# Returns the acceleration of the pitch angle of the inertial frame at:
# time "t" (in seconds)
# height "h" (in meters)
# velocity "v" (in m/s)
# pitch angle "pitch" (in radians measured from the vertical)
# change in pitch angle "v_pitch" (in radians)
# angle of the velocity vector "upsilon" (in radians measured from the vertical)
def eom_a_pitch(t, h, v, pitch, v_pitch, upsilon):
    Q = q(h, v)
    mach = get_mach(h, v)
    N = Q * S_area * aerodynamic_LUT(t, mach)[3] * (pitch - upsilon)
    a_pitch = (N * (aerodynamic_LUT(t, mach)[0]/100 - C_g_LUT(t)) + (Q * S_area * d**2 * pitch_damp * v_pitch)/(2 * v))/get_inertial_moment(t)
    return a_pitch

"""Functions that describe each RK45 step"""
# Returns the derivative of the state vector at the given time and state for the RK45 integrator
# time "t" (in seconds)
# state vector "x" containing [y, z, v_y, v_z, pitch, v_pitch, m]
def f(t, x):
    y, z, v_y, v_z, pitch, v_pitch, m = x

    v = np.sqrt(v_z**2 + v_y**2)

    a_y = eom_a_y(t, z, v, pitch, np.arctan2(v_y, v_z), m)
    a_z = eom_a_z(t, z, v, pitch, np.arctan2(v_y, v_z), m)
    a_pitch = eom_a_pitch(t, z, v, pitch, v_pitch, np.arctan2(v_y, v_z))

    return np.array([v_y, v_z, a_y, a_z, v_pitch, a_pitch, mass_flow_rate(t)])

# Returns the list of approximations at the given trial point for:
# the time at which the derivative is evaluated "trial_point" (in seconds)
# the step size for the finite difference approximation "step_size" (in seconds)
# the current state vector "state_vector" containing [y, z, v_y, v_z, pitch, v_pitch, m] in a numpy array
def rk_45_step(trial_point, step_size, state_vector):
    k1 = step_size * f(trial_point + step_size * b_tableau[0][0], state_vector)
    k2 = step_size * f(trial_point + step_size * b_tableau[1][0], state_vector + k1 * b_tableau[1][1])
    k3 = step_size * f(trial_point + step_size * b_tableau[2][0], state_vector + k1 * b_tableau[2][1] + k2 * b_tableau[2][2])
    k4 = step_size * f(trial_point + step_size * b_tableau[3][0], state_vector + k1 * b_tableau[3][1] + k2 * b_tableau[3][2] + k3 * b_tableau[3][3])
    k5 = step_size * f(trial_point + step_size * b_tableau[4][0], state_vector + k1 * b_tableau[4][1] + k2 * b_tableau[4][2] + k3 * b_tableau[4][3] + k4 * b_tableau[4][4])
    k6 = step_size * f(trial_point + step_size * b_tableau[5][0], state_vector + k1 * b_tableau[5][1] + k2 * b_tableau[5][2] + k3 * b_tableau[5][3] + k4 * b_tableau[5][4] + k5 * b_tableau[5][5])

    b4 = state_vector + k1 * b_tableau[6][1] + k2 * b_tableau[6][2] + k3 * b_tableau[6][3] + k4 * b_tableau[6][4] + k5 * b_tableau[6][5]
    b5 = state_vector + k1 * b_tableau[7][1] + k2 * b_tableau[7][2] + k3 * b_tableau[7][3] + k4 * b_tableau[7][4] + k5 * b_tableau[7][5] + k6 * b_tableau[7][6]
    return np.array(b4), np.array(b5)

# Returns a boolean and the normalized error denoting whether a step of the integrator failed or not due to significant normalized error for:
# some 4th order vector "b4"
# some 5th order vector "b5"
# the previous state "xi"
def check_step(b4, b5, xi):
    error_scale = atol_all + rtol * np.maximum(np.abs(xi), np.abs(b5))
    error = np.abs(b4 - b5)
    error_sum = 0
    N = 0
    while N < len(error):
        error_sum += (error[N]/error_scale[N]) ** 2
        N += 1
    normalized_error = np.sqrt(error_sum/N) # Normalizes the error between all values (regardless of units)
    if normalized_error <= 1:
        return True, normalized_error
    return False, normalized_error

# Returns the adjusted step size due to the amount of error in the previous step for:
# the previous step size "step_size"
# the normalized error "normalized_error"
def adjust_step(step_size, normalized_error):
    if normalized_error == 0:
        step_adjust = rk45_max_adjust # Increases step when there is no present error to save compute (less steps)
    else:
        step_adjust = rk45_step_constant * normalized_error ** (-1/5) # Modifies the step based on the current amount of normalized error
    step_adjust = min(rk45_max_adjust, max(rk45_min_adjust, step_adjust))
    return step_size * step_adjust # Add something in the integrator stack to check for too small or too large of step sizes
