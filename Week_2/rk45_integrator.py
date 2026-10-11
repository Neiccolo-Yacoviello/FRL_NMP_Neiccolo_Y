from integrator_functions import *

print(gravity(3000))

print(f"rho: {atmospheric_LUT(3000)[0]}, P: {atmospheric_LUT(3000)[1]}, T: {atmospheric_LUT(3000)[2]}")

print(weight(3000, m0))


