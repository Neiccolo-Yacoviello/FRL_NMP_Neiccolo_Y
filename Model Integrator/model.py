from model_functions import *

# Runs the integration mdoel with the PID controller for:
# The step size "step" in seconds
# The initial conditions "initial_condition" in the form of a numpy array
# The initial time "initial_time" in seconds
def model(step, initial_condition, initial_time, PID = True):
    if PID:
        reset_controller()
    state_vector = initial_condition
    previous_state = state_vector
    step_size = step
    time = initial_time
    altitude_list = [[0.0, 0.0]]
    v_list = [0.0]
    max_velo = 0
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
            v_list.append(velo)
            if PID:
                update_controller(time, state_vector, max_velo)

        max_velo = max(v_list)
        step_size = adjust_step(step_size, checked_step[1])
        if state_vector[3] > 0 and time > burnout_time:
            step_size = min(step_size, controller_step)

        safe = False # Recognizes whether or not the rocket hit the apogee target +/- 10m
        above = None # Recognizes the direction of the rocket overshoot (if applicable)
        if target_apogee + 10 > max(altitude_list)[0] > target_apogee - 10:
            safe = True
        if max(altitude_list)[0] < target_apogee - 10:
            above = False
        elif max(altitude_list)[0] > target_apogee + 10:
            above = True

        touchdown = False # Recognizes whether or not the rocket touched down
        if state_vector[1] < 0 and time > burnout_time:
            touchdown = True

        if safe and touchdown:
            # print(f"Simulation success; Rocket touch down at {time}s.")
            return max(altitude_list)[0], "success"
        
        elif np.abs(state_vector[4]) > np.radians(45) and 0.5 < time < 12.5:
            # print(f"Simulation ended; Rocket over-rotated at {time}s.")
            return max(altitude_list)[0], "over-rotation"
        
        elif above and touchdown:
            # print(f"Simulation ended; Rocket undershot the target apogee at {time}s.")
            return max(altitude_list)[0], "overshoot"

        elif not above and touchdown:
            # print(f"Simulation ended; Rocket overshot the target apogee at {time}s.")
            return max(altitude_list)[0], "undershoot"