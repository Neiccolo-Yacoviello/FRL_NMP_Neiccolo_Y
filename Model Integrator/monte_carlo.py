from model import *
import random

gains = [1.0, 0.0, 0.2] # The gains for the PID controller (K_P, K_I, K_D)

seed = 1 # The seed for the simulation analysis, used to ensure reproducibility of the simulation results
uncertainty_limit = 0.025 # The limit of uncertainty for inital conditions (such as wind, mass, and thrust factors), +/- x.x%
wind_limit = 5 # The limit of wind speed in m/s, +/- x.x m/s
trials = 100 # The number of trials to run for the monte carlo analysis

initial_conditions = np.array([0, 0, 0, 0, 0, 0, m0 * MCV["mass_offset"]])

random.seed(seed)

sim_success_list = []
sim_fail_list = []

over_rotates = 0
undershoots = 0
overshoots = 0

x = 0
while x < trials:
    wind_speed = round(random.uniform(-wind_limit, wind_limit), 2) # The wind speed in m/s
    wind_direction = round(random.uniform(-180, 180), 2) # The wind direction in degrees
    thrust_scale = round(random.uniform(1 - uncertainty_limit, 1 + uncertainty_limit), 4) # The scale factor for the thrust of the rocket engine
    burn_scale = round(random.uniform(1 - uncertainty_limit, 1 + uncertainty_limit), 4) # The scale factor for the burn time of the rocket engine
    drag_scale = round(random.uniform(1 - uncertainty_limit, 1 + uncertainty_limit), 4) # The scale factor for the drag of the rocket
    mass_offset = round(random.uniform(1 - uncertainty_limit, 1 + uncertainty_limit), 4) # The scale factor for the mass of the rocket

    MCV["wind_speed"] = wind_speed
    MCV["wind_direction"] = wind_direction
    MCV["thrust_scale"] = thrust_scale
    MCV["burn_scale"] = burn_scale
    MCV["drag_scale"] = drag_scale
    MCV["mass_offset"] = mass_offset

    max_info = model(0.01, initial_conditions, 0)
    print([round(max_info[0], 3), max_info[1]])

    if max_info[1] == "success":
        sim_success_list.append(round(max_info[0], 3))
    elif max_info[1] == "over-rotation":
       sim_fail_list.append((round(max_info[0], 3), [round(val, 3) for val in list(MCV.values())], max_info[1]))
       over_rotates += 1
    elif max_info[1] == "undershoot":
        sim_fail_list.append((round(max_info[0], 3), [round(val, 3) for val in list(MCV.values())], max_info[1]))
        undershoots += 1
    elif max_info[1] == "overshoot":
        sim_fail_list.append((round(max_info[0], 3), [round(val, 3) for val in list(MCV.values())], max_info[1]))
        overshoots += 1

    x += 1

print(f"\n\nNumber of failed simulations: {len(sim_fail_list)}")
print(f"Number of successful simulations: {len(sim_success_list)}")
print(f"Max apogee of failed simulations: {max(sim_fail_list, key=lambda x: x[0])[0]}m")
print(f"Min apogee of failed simulations: {min(sim_fail_list, key=lambda x: x[0])[0]}m")
print(f"Number of over-rotation failures: {over_rotates}")
print(f"Number of undershoot failures: {undershoots}")
print(f"Number of overshoot failures: {overshoots}")

for i in sim_fail_list:
    print(f"Failed simulation apogee: {i[0]}m, with initial conditions: {i[1]}, in this direction: {i[2]}")