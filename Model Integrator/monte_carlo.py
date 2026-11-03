from model import *
import random

gains = [3.0, 0.25, 0.95] # The gains for the PID controller (K_P, K_I, K_D)

seed = 5421906848956 # The seed for the simulation analysis, used to ensure reproducibility of the simulation results
thrust_uncertainty_limit = 0.05 # The limit of uncertainty for thrust (I_sp) +/- x.x%
burn_uncertainty_limit = 0.05 # The limit of uncertainty for burn time (burnout_time) +/- x.x%
drag_uncertainty_limit = 0.04 # The limit of uncertainty for drag (C_d) +/- x.x%
mass_uncertainty_limit = 0.1 # The limit of uncertainty for mass (m) +/- x.x%
wind_limit = 5 # The limit of wind speed +/- x.x m/s
trials = 2000 # The number of trials to run for the monte carlo analysis

initial_conditions = np.array([0, 0, 0, 0, 0, 0, m0 * MCV["mass_offset"]])

random.seed(seed)

sim_success_list = []
sim_fail_list = []

over_rotates = 0
undershoots = 0
failures = 0
overshoots = 0

overshoot_conditions_avg = [0, 0, 0, 0, 0, 0]
undershoot_conditions_avg = [0, 0, 0, 0, 0, 0]

x = 0
while x < trials:
    wind_speed = round(random.uniform(-wind_limit, wind_limit), 2) # The wind speed in m/s
    wind_direction = round(random.uniform(-180, 180), 2) # The wind direction in degrees
    thrust_scale = round(random.uniform(1 - thrust_uncertainty_limit, 1 + thrust_uncertainty_limit), 5) # The scale factor for the thrust of the rocket engine
    burn_scale = round(random.uniform(1 - burn_uncertainty_limit, 1 + burn_uncertainty_limit), 5) # The scale factor for the burn time of the rocket engine
    drag_scale = round(random.uniform(1 - drag_uncertainty_limit, 1 + drag_uncertainty_limit), 5) # The scale factor for the drag of the rocket
    mass_offset = round(random.uniform(1 - mass_uncertainty_limit, 1 + mass_uncertainty_limit), 5) # The scale factor for the mass of the rocket

    MCV["wind_speed"] = wind_speed
    MCV["wind_direction"] = wind_direction
    MCV["thrust_scale"] = thrust_scale
    MCV["burn_scale"] = burn_scale
    MCV["drag_scale"] = drag_scale
    MCV["mass_offset"] = mass_offset

    x += 1

    max_info = model(0.01, initial_conditions, 0)
    # print([float(round(max_info[0], 3)), max_info[1], x])

    if max_info[1] == "success":
        sim_success_list.append(round(max_info[0], 3))
    elif max_info[1] == "over-rotation":
       sim_fail_list.append((round(max_info[0], 3), [val for val in list(MCV.values())], max_info[1], x))
       over_rotates += 1
    elif -1 < max_info[0] < 1:
        print(f" Rocket failed to launch, with initial conditions: {list(MCV.values())}")
        failures += 1
    elif max_info[1] == "undershoot":
        sim_fail_list.append((round(max_info[0], 3), [val for val in list(MCV.values())], max_info[1], x))
        undershoots += 1
    elif max_info[1] == "overshoot":
        sim_fail_list.append((round(max_info[0], 3), [val for val in list(MCV.values())], max_info[1], x))
        overshoots += 1

print(f"\n\nNumber of failed simulations: {len(sim_fail_list)}")
print(f"Number of successful simulations: {len(sim_success_list)}, success rate: {round(len(sim_success_list)/trials * 100, 2)}%")
print(f"Average apogee of successful simulations: {round(sum(sim_success_list)/len(sim_success_list), 3)}m")
print(f"\nMax apogee of failed simulations: {max(sim_fail_list, key=lambda x: x[0])[0]}m")
print(f"Min apogee of failed simulations: {min(sim_fail_list, key=lambda x: x[0])[0]}m")
print(f"\nNumber of over-rotation failures: {over_rotates}")
print(f"Number of undershoot failures: {undershoots}")
print(f"Number of overshoot failures: {overshoots}")

for i in sim_fail_list:
    # print(f"Failed simulation apogee: {i[0]}m, with initial conditions: {i[1]}, in this direction: {i[2]}, for trial number: {i[3]}")
    if i[2] == "overshoot":
        for j in i[1]:
            overshoot_conditions_avg[i[1].index(j)] += j
    elif i[2] == "undershoot":
        for j in i[1]:
            undershoot_conditions_avg[i[1].index(j)] += j

print("\n[thrust_scale, burn_scale, drag_scale, mass_offset]")
if overshoots > 0:
    print(f"Average conditions for overshoot failures: { [round(i/overshoots, 5) for i in overshoot_conditions_avg[2:]] }")
if undershoots > 0:
    print(f"Average conditions for undershoot failures: { [round(i/undershoots, 5) for i in undershoot_conditions_avg[2:]] }")