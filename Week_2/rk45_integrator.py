from integrator_functions import *
from typing import Final

print(f"Gravity at h=3000: {gravity(3000)}")

print(f"Density at h=3000: {atmospheric_LUT(3000)[0]}")
print(f"Pressure at h=3000: {atmospheric_LUT(3000)[1]}")
print(f"Temperature at h=3000: {atmospheric_LUT(3000)[2]}")

print(f"Weight at h=3000: {weight(3000, 4.1)}")

print(f"Thrust at t=3.0: {thrust_LUT(3.0)}")

print(f"Mass flow at t=2.9: {mass_flow_rate(2.9)}")

print(f"Center of gravity at t=2.0: {C_g_LUT(2.0)}")

print(f"Dynamic Pressure (Q) at h=1000, v=250: {q(1000, 250)}")

print(f"Aerodynamic info [Z_cp, C_d, C_A, C_N]: {aerodynamic_LUT(3, get_mach(1000, 250))}")

print(f"Coefficient of Lift: {get_C_L(3, 1000, 250, 0.1, 0.2)}")

print(f"Lift: {lift(3, 1000, 250, 0.1, 0.2)}")

print(f"Drag: {drag(3, 1000, 250)}")

print(f"Acceleration in \"z\": {eom_a_z(3, 1000, 250, 0.1, 0.2, 4.1)}")

print(f"Acceleration in \"y\": {eom_a_y(3, 1000, 250, 0.1, 0.2, 4.1)}")

print(f"Acceleration in \"pitch\": {eom_a_pitch(3, 1000, 250, 0.1, 0.01, 0.2)}")
