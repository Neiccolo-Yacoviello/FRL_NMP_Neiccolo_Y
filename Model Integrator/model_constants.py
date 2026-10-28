from typing import Final
from pathlib import Path
import numpy as np

DATA_DIR = Path(__file__).parent

"""APOGEE TARGET"""
target_apogee = 3000

"""EOM Constants"""
g_Earth_S: Final = 9.81 # The acceleration due to gravity at Earth's surface

R_Earth: Final = 6371008 # The mean radius of the Earth in meters

I_sp: Final = 2284 # The specific impulse of the rocket engine in Newton seconds

air_specific_heat_ratio: Final = 1.4 # The standard specific heat ratio for dry air

air_specific_gas_constant: Final = 287.05 # The specific gas constant for dry air

m0: Final = 4.630 # The initial mass of the rocket in kilograms

m_dry: Final = 3.338 # The mass of the rocket after its burn in kilograms

d: Final = 0.0762 # The diameter of the rocket in meters

L: Final = 1.61 # The length of the rocket in meters

S_area: Final = np.pi * (d / 2)**2 # The reference area of the rocket in m^2

pitch_damp: Final = -0.0005 # The pitch dampening coefficient

burnout_time: Final = 3.593 # The time at which the motor burns out

C_Na: Final = 7.392 # The Normal Force per radian of the rocket

max_velo: Final = 417.5 # The max speed in m/s the rocket reaches during its flight (in perfect conditions)

drogue_chute_d: Final = 0.457 # The diameter of the drogue chute in meters

main_chute_d: Final = 1.0 # The diameter of the main chute in meters

drogue_area: Final = np.pi * (drogue_chute_d / 2)**2 # The area of the drogue chute in m^2

main_area: Final = np.pi * (main_chute_d / 2)**2 # The area of the main chute in m^2

chute_C_d: Final = 0.8 # The drag coefficient of the parachutes

ab_L: Final = 0.0493 # The length of the air brake panels

ab_W: Final = 0.05388 # The width of the air brake panels

ab_C_d: Final = 1.3 # The drag coefficient of the air brake panels

"""Monte Carlo Analysis and PID constants"""
MCV = { # The variables used during monte carlo analysis to slightly vary initial conditions
    "wind_speed": 0.0,
    "wind_direction": 0.0,
    "thrust_scale": 1.0,
    "burn_scale": 1.0,
    "drag_scale": 1.0,
    "mass_offset": 0.0,
}

gains = [0.5, 1, 0] # The gains for the PID controller (K_P, K_I, K_D)

"""Raw CSV data"""
thrust_data: Final = np.loadtxt(DATA_DIR / "thrust(updated_10_26).csv", delimiter=',', skiprows=1) # The raw thrust data

c_g_data: Final = np.loadtxt(DATA_DIR / "c_g(updated_10_26).csv", delimiter=',', skiprows=1) # The raw cg data

raw_A_data = np.loadtxt(DATA_DIR / "aerodynamics(updated_10_26).csv", delimiter=',', skiprows=1) # The raw aerodynamics data
burnout = raw_A_data[:, 3].argmax() # Largest mach value
chute = np.argmax(raw_A_data[:, 1] > 10) # Jump in the drag coefficient
aerodynamics_data: Final = { # Describes the main phases of flight to pull aero data from since "every" mach number is reached twice on ascent (3 times for smaller values with descent)
    "boost": raw_A_data[:burnout + 1][np.argsort(raw_A_data[:burnout + 1, 3])],
    "coast": raw_A_data[burnout:chute][np.argsort(raw_A_data[burnout:chute, 3])],
    "descent": raw_A_data[chute:][np.argsort(raw_A_data[chute:, 3])],
}

"""EOM Look Up Tables (from CSV data)"""
# Returns the atmospheric properties (density, pressure, and temperature) at:
# height "h" (in meters) above the Earth's surface
def atmospheric_LUT(h):
    if h >= 11000 and h < 25000:
        T = -56.46
        P = 22.65 * np.exp(1.73 - 0.000157 * h)
    else:
        T = 15.04 - 0.00649 * h
        P = 101.29 * ((T + 273.1) / 288.08)**5.256
    rho = P / (0.2869 * (T + 273.1))    
    return rho, P, T

# Returns the aerodynamic info is given by [Z_cp, C_d, C_A, C_N, v]
# some velocity "v" (in mach number)
# some velocity in the z direction "v_z" (in m/s)
def aerodynamic_LUT(t, v, v_z):
    phase = "boost"
    if t > burnout_time: 
        phase = "coast"
        if v_z < 0:
            phase = "descent"
    data = aerodynamics_data[phase]
    v_data = data[:, 3]
    return [float(np.interp(v, v_data, data[:, i])) for i in range(data.shape[1] - 1)]

# Returns the interpolated thrust of the engine at:
# some time "t" (in seconds) where after burnout (t > 3.593 seconds) the thrust is zero
def thrust_LUT(t):
    t = t / MCV["burn_scale"]
    if t < burnout_time: # Burnout of the current motor (K700)
        return np.interp(t, thrust_data[:, 0], thrust_data[:, 1]) * MCV["thrust_scale"]
    else:
        return 0

# Returns the interpolated center of gravity of the rocket (from the nose) at:
# some time "t" (in seconds) where after burnout (t > 3.593 seconds) the center of gravity is constant
def C_g_LUT(t):
    if t < burnout_time: # Burnout of the current motor (K700)
        return np.interp(t, c_g_data[:, 0], c_g_data[:, 1])
    else:
        return c_g_data[-1, 1]  # Return the last known center of gravity

"""RK45 Constants"""
b_tableau: Final = [
    [0.0  , 0.0      , 0.0       , 0.0       , 0.0        , 0.0   , 0.0 ],
    [0.25 , 0.25     , 0.0       , 0.0       , 0.0        , 0.0   , 0.0 ],
    [3/8  , 3/32     , 9/32      , 0.0       , 0.0        , 0.0   , 0.0 ],
    [12/13, 1932/2197, -7200/2197, 7296/2197 , 0.0        , 0.0   , 0.0 ],
    [1.0  , 439/216  , -8.0      , 3680/513  , -845/4104  , 0.0   , 0.0 ],
    [0.5  , -8/27    , 2.0       , -3544/2565, 1859/4104  , -11/40, 0.0 ],
    [0.0  , 25/216   , 0.0       , 1408/2565 , 2197/4104  , -1/5  , 0.0 ], # B4 coefficients for RK45 (4th order)
    [0.0  , 16/135   , 0.0       , 6656/12825, 28561/56430, -9/50 , 2/55], # B5 coefficients for RK45 (5th order)
]

rk45_step_constant: Final = 0.84 # The step size adjustment factor for when a potential rk45 solution is rejected due to truncation error

rk45_max_adjust: Final = 2 # The max step size scalar allowed

rk45_min_adjust: Final = 0.5 # The min step size scalar allowed

atol: Final = 1e-4 # The allowable absolute error for any value

atol_all: Final = (1e-3, 1e-3, 1e-6, 1e-6, 1e-3, 1e-3, 1e-6) # The allowable absolute error for each value in the state vector

rtol: Final = 1e-3 # The allowable relative error for any value

"""Moments of inertia in each axes of rotation (only x included)"""
def get_inertial_moment(t): # Changes during the burn
    inertial_moment_x  = 0.4422 # Max exported from OR
    if t < burnout_time:
        return inertial_moment_x + (0.1780 * (1 - t / burnout_time)) # .1780 is the difference between the max and min exported from OR
    return inertial_moment_x