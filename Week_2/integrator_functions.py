from integrator_constants import *

# Returns the acceleration due to gravity given some height "h" above the Earth's surface (for it's mean radius)
def get_gravity(h):
    g = g_Earth_S * (R_Earth / (R_Earth + h))**2
    return g

