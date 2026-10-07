from typing import Final

"""EOM Constants"""
g_Earth_S: Final = 9.81 # The acceleration due to gravity at Earth's surface

R_Earth: Final = 6371008 # The mean radius of the Earth in meters

I_sp: Final = 300 # The specific impulse of the rocket engine in seconds

m0: Final = 5.5 # The initial mass of the rocket in kilograms

"""EOM Look Up Tables"""
atmospheric_LUT: Final = [ # For some height "h" (in meters), the atmospheric info is given by [h, rho, P, T]

]

aerodynamic_LUT: Final = [ # For some velocity "v" (in m/s / mach number), the aerodynamic info is given by [v, C_d, C_N, Z_cp]

]

thrust_LUT: Final = [ # For some time "t" (in seconds), the thrust info is given by [t, T]

]

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

truncation_error: Final = 1e-6 # The lower limit for accepable truncation error in the rk45 method

"""Moments of inertia in each axes of rotation"""
inertia_tensor: Final = [
    ["""I_xx""", 0.0, 0.0],
    [0.0, """I_yy""", 0.0],
    [0.0, 0.0, """I_zz"""]
]
