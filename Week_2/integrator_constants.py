from typing import Final

g_Earth_S: Final = 9.81 # the acceleration due to gravity at Earth's surface

R_Earth: Final = 6371008 # the mean radius of the Earth in meters

b_tableau: Final = [
    [0.0  , 0.0      , 0.0       , 0.0       , 0.0        , 0.0   , 0.0 ],
    [0.25 , 0.25     , 0.0       , 0.0       , 0.0        , 0.0   , 0.0 ],
    [3/8  , 3/32     , 9/32      , 0.0       , 0.0        , 0.0   , 0.0 ],
    [12/13, 1932/2197, -7200/2197, 7296/2197 , 0.0        , 0.0   , 0.0 ],
    [1.0  , 439/216  , -8.0      , 3680/513  , -845/4104  , 0.0   , 0.0 ],
    [0.5  , -8/27    , 2.0       , -3544/2565, 1859/4104  , -11/40, 0.0 ],
    [0.0  , 16/135   , 0.0       , 6656/12825, 28561/56430, -9/50 , 2/55], # B4 coefficients for RK45 (4th order)
    [0.0  , 25/216   , 0.0       , 1408/2565 , 2197/4104  , -1/5  , 0.0 ], # B5 coefficients for RK45 (5th order)
]









truncation_error: Final = 1e-6

