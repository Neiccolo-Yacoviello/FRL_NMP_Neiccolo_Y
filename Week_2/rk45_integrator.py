from integrator_functions import *
from integrator_constants import *

initial_conditions = np.array([0, 0, 0, 0, 0, 0, m0]) # the beginning state vector containing [y, z, v_y, v_z, pitch, v_pitch, m]

def model(step, initial):
    state_vector = initial
    previous_state = initial
    step_size = step
    time = 0
    while True:
        diff = rk_45_step(time, step_size, state_vector) # Contains the b4 and b5 evaluations of the step
        checked_step = check_step(diff[0], diff[1], previous_state) # Contains the step verification and the error of the step
        if checked_step[0]:
            previous_state = state_vector
            state_vector = diff[1]
            print(f"Z: {round(state_vector[1], 3)}, V_Z: {round(state_vector[3], 3)}, V_Y: {round(state_vector[2], 3)}, Pitch: {round(state_vector[4], 3)}, Mass: {round(state_vector[6], 3)}, t: {round(time, 4)}")
            time += step_size

        step_size = adjust_step(step_size, checked_step[1])

        if state_vector[1] < 0:
            print("Simulation ended; Rocket touched down.")
            break
        elif np.abs(state_vector[4]) > 30 and time < 25: # Temporary (might remove later)
            print("Simulation ended; Rocket over-rotated.")
            break

model(0.01, initial_conditions)