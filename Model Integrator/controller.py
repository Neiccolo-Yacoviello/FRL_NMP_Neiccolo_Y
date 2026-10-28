from model_constants import *
import numpy as np

controller_step = 0.05

controller_state = {
    "control_value": 0.0,
    "next_step": 0.0,
    "integral": 0.0,
    "previous_error": None,
    "predicted_apogee": None,
}

# Uses a simple Eulers method to predict the apogee for the PID controller for:
# The current time of flight "t" (in s)
# The current height of the rocket  "h" (in m)
# The current velocity in the y direction "v_y" (in m/s)
# The current velocity in the z direction "v_z" (in m/s)
# The current mass of the rocket "m" (in kg; constant after burnout/always constant in this function)
# The value of the PID controller for the air brakes "control value" (0 <= value <= 1)
# The step size "step" (in s)
def predict_apogee(t, h, v_y, v_z, m, control_value, step):
    while v_z > 0:
        v = np.sqrt(v_y**2 + v_z**2)
        temp = atmospheric_LUT(h)[2] + 273.15
        rho = atmospheric_LUT(h)[0]
        local_mach = np.sqrt(temp * air_specific_gas_constant * air_specific_heat_ratio)
        mach = v/local_mach
        C_d = aerodynamic_LUT(t, mach, v_z)[1]
        ab_area = 4 * ab_W * ab_L * np.sin(np.radians(20 * control_value))

        D = 0.5 * rho * v**2 * (C_d * S_area + ab_C_d * ab_area)

        g = g_Earth_S * (R_Earth / (R_Earth + h))**2

        v_y -= step * (D * v_y/v)/m  
        v_z -= step * ((D * v_z/v)/m + g)

        h += v_z * step

    return h

# Returns the signal of the PID for:
# The target apogee "set_point" (3000m)
# The height of the rocket at some point time "value"
# The step size taken during the process "step" (in s)
def pid(set_point, value, step):
    K_P, K_I, K_D = gains
    error = value - set_point
    if controller_state["previous_error"] is None:
        derivative = 0.0
    else:
        derivative = (error - controller_state["previous_error"])/step
    controller_state["previous_error"] = error
    signal = K_P * error + K_I * (controller_state["integral"] + error * step) + K_D * derivative

    if 0.0 < signal and signal < 1.0: 
        controller_state["integral"] += error * step

    return min(1.0, max(signal, 0.0))

# Updates the apogee prediction and PID signal for the drag function to read for EOMs
def update_controller(time, state):
    if time < controller_state["next_step"]:
        return
    controller_state["next_step"] = time + controller_step
    
    y, z, v_y, v_z, pitch, v_pitch, m = state
    v = np.sqrt(v_y**2 + v_z**2)

    if time > burnout_time and v_z > 0 and 0.25 * max_velo < v and v < 0.5 * max_velo:
        controller_state["predicted_apogee"] = predict_apogee(time, z, v_y, v_z, m, 0, controller_step)
        controller_state["control_value"] = pid(target_apogee, controller_state["predicted_apogee"], controller_step)
    else:
        controller_state["control_value"] = 0.0

# Resets the controller for another run of the model
def reset_controller():
    controller_state.update(control_value = 0.0, next_step = 0.0, integral = 0.0, previous_error = None, predicted_apogee = None)