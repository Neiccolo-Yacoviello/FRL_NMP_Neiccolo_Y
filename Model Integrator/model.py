from model_functions import *

initial_conditions = np.array([0, 0, 0, 0, 0, 0, m0 + MCV["mass_offset"]]) # the beginning state vector containing [y, z, v_y, v_z, pitch, v_pitch, m]

def model(step, initial_condition, initial_time):
    reset_controller()
    state_vector = initial_condition
    previous_state = state_vector
    step_size = step
    time = initial_time
    altitude_list = []
    while True:
        diff = rk_45_step(time, step_size, state_vector) # Contains the b4 and b5 evaluations of the step
        checked_step = check_step(diff[0], diff[1], previous_state) # Contains the step verification and the error of the step
        velo = np.sqrt(state_vector[2]**2 + state_vector[3]**2)
        if checked_step[0]:
            previous_state = state_vector
            state_vector = diff[1]
            time += step_size
            # print(f"Z: {round(state_vector[1], 3)}, V_Z: {round(state_vector[3], 3)}, V_Y: {round(state_vector[2], 3)}, Pitch: {round(state_vector[4], 3)}, Mass: {round(state_vector[6], 3)}, t: {round(time, 4)}")
            altitude_list.append([state_vector[1], time])
            update_controller(time, state_vector)

        step_size = adjust_step(step_size, checked_step[1])
        if state_vector[3] > 0:
            step_size = min(step_size, controller_step)
     
        if state_vector[1] < 0:
            # print(f"Simulation ended; Rocket touched down at {time}s.")
            return max(altitude_list)[0], max(altitude_list)[1], time, True
        elif np.abs(state_vector[4]) > 1 and time < burnout_time:
            # print(f"Simulation ended; Rocket over-rotated at {time}s.")
            return max(altitude_list)[0], max(altitude_list)[1], time, False



max_info = model(0.01, initial_conditions, 0)

print(f"The rocket reached apogee at {round(max_info[0], 3)}m at {round(max_info[1], 3)}s and touched down at {round(max_info[2], 3)}s")