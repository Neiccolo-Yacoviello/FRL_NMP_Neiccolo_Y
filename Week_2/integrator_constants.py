from typing import Final
from pathlib import Path
import numpy as np

DATA_DIR = Path(__file__).parent

"""EOM Constants"""
g_Earth_S: Final = 9.81 # The acceleration due to gravity at Earth's surface

R_Earth: Final = 6371008 # The mean radius of the Earth in meters

I_sp: Final = 2284 # The specific impulse of the rocket engine in Newton seconds

m0: Final = 4.637 # The initial mass of the rocket in kilograms

d: Final = 0.0762 # The diameter of the rocket in meters

L: Final = 1.56 # The length of the rocket in meters

nose_C_p: Final = 1.12 # The center of pressure measued from the nose of the rocket in meters

S_area: Final = np.pi * (d / 2)**2 # The reference area of the rocket in m^2

C_m_q: Final = -0.5 # The pitch dampening coefficient

thrust_data: Final = np.loadtxt(DATA_DIR / "thrust.csv", delimiter=',', skiprows=1) # The raw thrust data

c_g_data: Final = np.loadtxt(DATA_DIR / "c_g.csv", delimiter=',', skiprows=1) # The raw cg data


"""EOM Look Up Tables"""
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

aerodynamic_LUT: Final = [ # For some velocity "v" (in m/s / mach number), the aerodynamic info is given by [v, C_d, C_L, C_N, Z_cp]

    # return np.array()
]

# Returns the interpolated thrust of the engine at:
# some time "t" (in seconds) where after burnout (t > 3.593 seconds) the thrust is zero
def thrust_LUT(t):
    if t < 3.593: # Burnout of the current motot (K700)
        return np.interp(t, thrust_data[:, 0], thrust_data[:, 1])
    else:
        return 0

# Returns the interpolated center of gravity of the rocket (from the nose) at:
# some time "t" (in seconds) where after burnout (t > 3.593 seconds) the center of gravity is constant
def C_g_LUT(t):
    if t < 3.593: # Burnout of the current motot (K700)
        return np.interp(t, c_g_data[:, 0], c_g_data[:, 1])/100
    else:
        return c_g_data[-1, 1]/100  # Return the last known center of gravity

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

truncation_error_upper: Final = 1e-3 # The upper limit for accepable truncation error in the rk45 method

"""Moments of inertia in each axes of rotation"""
inertia_tensor: Final = [
    ["""I_xx""", 0.0, 0.0],
    [0.0, """I_yy""", 0.0],
    [0.0, 0.0, """I_zz"""]
]
