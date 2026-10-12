from integrator_functions import *

print(f"Gravity at h=3000: {gravity(3000)}")

print(f"Density at h=3000: {atmospheric_LUT(3000)[0]}")
print(f"Pressure at h=3000: {atmospheric_LUT(3000)[1]}")
print(f"Temperature at h=3000: {atmospheric_LUT(3000)[2]}")

print(f"Weight at h=3000: {weight(3000, m0)}")

print(f"Thrust at t=3.0: {thrust_LUT(3.0)}")

print(f"Center of gravity at t=2.0: {C_g_LUT(2.0)}")
