from model import *
import random
import time
from multiprocessing import Pool

gains_monte_carlo = [2.07, 0.0, 0.76] # The gains for the PID controller (K_P, K_I, K_D)

seed = 3247964780731 # The seed for the simulation analysis, used to ensure reproducibility of the simulation results
thrust_uncertainty_limit = 0.05 # The limit of uncertainty for thrust (I_sp) +/- x.x%
burn_uncertainty_limit = 0.05 # The limit of uncertainty for burn time (burnout_time) +/- x.x%
drag_uncertainty_limit = 0.04 # The limit of uncertainty for drag (C_d) +/- x.x%
mass_uncertainty_limit = 0.1 # The limit of uncertainty for mass (m) +/- x.x%
wind_limit = 5 # The limit of wind speed +/- x.x m/s
trials = 10000 # The number of trials to run for the monte carlo analysis

sim_success_list = []
sim_fail_list = []

over_rotates = 0
undershoots = 0
failures = 0
overshoots = 0

overshoot_conditions_avg = [0, 0, 0, 0, 0, 0]
undershoot_conditions_avg = [0, 0, 0, 0, 0, 0]

def run_trial(job):
    gains[:] = gains_monte_carlo
    z, parameters = job
    MCV.update(parameters)
    initial_conditions = np.array([0, 0, 0, 0, 0, 0, m0 * MCV["mass_offset"]])
    max_info = model(0.01, initial_conditions, 0)
    return z, max_info, list(MCV.values())

if __name__ == "__main__":
    start = time.perf_counter()
    random.seed(seed)

    jobs = []

    for i in range(1, trials + 1):
        parameters = {
        "wind_speed": round(random.uniform(-wind_limit, wind_limit), 2), # The wind speed in m/s
        "wind_direction": round(random.uniform(-180, 180), 2), # The wind direction in degrees
        "thrust_scale": round(random.uniform(1 - thrust_uncertainty_limit, 1 + thrust_uncertainty_limit), 5), # The scale factor for the thrust of the rocket engine
        "burn_scale": round(random.uniform(1 - burn_uncertainty_limit, 1 + burn_uncertainty_limit), 5), # The scale factor for the burn time of the rocket engine
        "drag_scale": round(random.uniform(1 - drag_uncertainty_limit, 1 + drag_uncertainty_limit), 5), # The scale factor for the drag of the rocket
        "mass_offset": round(random.uniform(1 - mass_uncertainty_limit, 1 + mass_uncertainty_limit), 5), # The scale factor for the mass of the rocket
        }
        jobs.append((i, parameters))


    with Pool() as pool:
        for z, max_info, conditions in pool.imap(run_trial, jobs):
            # print([float(round(max_info[0], 3)), max_info[1], z])

            if max_info[1] == "success":
                sim_success_list.append(round(max_info[0], 3))
            elif max_info[1] == "over-rotation":
                sim_fail_list.append((round(max_info[0], 3), conditions, max_info[1], z))
                over_rotates += 1
            elif -1 < max_info[0] < 1:
                print(f" Rocket failed to launch, with initial conditions: {conditions}")
                failures += 1
            elif max_info[1] == "undershoot":
                sim_fail_list.append((round(max_info[0], 3), conditions, max_info[1], z))
                undershoots += 1
            elif max_info[1] == "overshoot":
                sim_fail_list.append((round(max_info[0], 3), conditions, max_info[1], z))
                overshoots += 1

    end = time.perf_counter()

    print(f"\n\nProgram took {(end-start):.2f}s to execute")
    print(f"K_P: {gains_monte_carlo[0]}, K_I: {gains_monte_carlo[1]}, K_D: {gains_monte_carlo[2]}")
    print(f"Number of successful simulations: {len(sim_success_list)}, success rate: {round(len(sim_success_list)/trials * 100, 2)}%")
    print(f"Average apogee of successful simulations: {round(sum(sim_success_list)/len(sim_success_list), 3)}m")
    # print(f"\nMax apogee of failed simulations: {max(sim_fail_list, key=lambda x: x[0])[0]}m")
    # print(f"Min apogee of failed simulations: {min(sim_fail_list, key=lambda x: x[0])[0]}m")
    print(f"\nNumber of over-rotation failures: {over_rotates}")
    print(f"Number of undershoot failures: {undershoots}")
    print(f"Number of overshoot failures: {overshoots}")

    for x in sim_fail_list:
        # print(f"Failed simulation apogee: {i[0]}m, with initial conditions: {i[1]}, in this direction: {i[2]}, for trial number: {i[3]}")
        if x[2] == "overshoot":
            for y in x[1]:
                overshoot_conditions_avg[x[1].index(y)] += y
        elif x[2] == "undershoot":
            for y in x[1]:
                undershoot_conditions_avg[x[1].index(y)] += y

    print("\n[thrust_scale, burn_scale, drag_scale, mass_offset]")
    if overshoots > 0:
        print(f"Average conditions for overshoot failures: {[round(x/overshoots, 3) for x in overshoot_conditions_avg[2:]]}")
    if undershoots > 0:
        print(f"Average conditions for undershoot failures: {[round(x/undershoots, 3) for x in undershoot_conditions_avg[2:]]}")