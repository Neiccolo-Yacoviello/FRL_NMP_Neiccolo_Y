from integrator_functions import *
from integrator_constants import *

initial_conditions = np.array([0, 0, 0, 0, 0, 0, m0]) # the beginning state vector containing [y, z, v_y, v_z, pitch, v_pitch, m]

def model(step, initial_condition, initial_time):
    state_vector = initial_condition
    previous_state = state_vector
    step_size = step
    time = initial_time
    velo_list = []
    altitude_list = []
    while True:
        diff = rk_45_step(time, step_size, state_vector) # Contains the b4 and b5 evaluations of the step
        checked_step = check_step(diff[0], diff[1], previous_state) # Contains the step verification and the error of the step
        if checked_step[0]:
            previous_state = state_vector
            state_vector = diff[1]
            # print(f"Z: {round(state_vector[1], 3)}, V_Z: {round(state_vector[3], 3)}, V_Y: {round(state_vector[2], 3)}, Pitch: {round(state_vector[4], 3)}, Mass: {round(state_vector[6], 3)}, t: {round(time, 4)}")
            velo_list.append([np.sqrt(state_vector[2]**2 + state_vector[3]**2), time])
            altitude_list.append([state_vector[1], time])
            time += step_size

        step_size = adjust_step(step_size, checked_step[1])

        # State Machine; use TBD
        # if time == 0: 
        #     state = "pad"
        # elif time < burnout_time:
        #     state = "boost"
        # elif time > burnout_time:
        #     state = "coast"
        # elif state_vector[3] < 0:
        #     state = "descent"

        if state_vector[1] < 0:
            print(f"Simulation ended; Rocket touched down at {time}s.")
            return max(velo_list), max(altitude_list)
        elif np.abs(state_vector[4]) > 1 and time < burnout_time:
            print(f"Simulation ended; Rocket over-rotated at {time}s.")
            break

max_info = model(0.1, initial_conditions, 0)
print(f"Max Velocity: {max_info[0][0]} at {max_info[0][1]}s")
print(f"Max Altitude: {max_info[1][0]} at {max_info[1][1]}s")