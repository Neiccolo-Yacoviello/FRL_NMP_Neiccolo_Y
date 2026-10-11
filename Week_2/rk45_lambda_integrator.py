from integrator_constants import *

# Returns the value of the derivative of the input function at:
# initial position "x"
# initial position "y"
# somve given lambda function f
def func(x, y, f):
    y_eval = f(x, y)
    return y_eval

# Returns the list of derivative approximations at:
# the point that the derivative is being evaluated at "trial_point" (x, y)
# the step size for the finite difference approximation "step_size" 
def rk_45_step(x, y, step_size, f):
    k1 = step_size * func(x + step_size * b_tableau[0][0], y, f)
    k2 = step_size * func(x + step_size * b_tableau[1][0], y + k1 * b_tableau[1][1], f)
    k3 = step_size * func(x + step_size * b_tableau[2][0], y + k1 * b_tableau[2][1] + k2 * b_tableau[2][2], f)
    k4 = step_size * func(x + step_size * b_tableau[3][0], y + k1 * b_tableau[3][1] + k2 * b_tableau[3][2] + k3 * b_tableau[3][3], f)
    k5 = step_size * func(x + step_size * b_tableau[4][0], y + k1 * b_tableau[4][1] + k2 * b_tableau[4][2] + k3 * b_tableau[4][3] + k4 * b_tableau[4][4], f)
    k6 = step_size * func(x + step_size * b_tableau[5][0], y + k1 * b_tableau[5][1] + k2 * b_tableau[5][2] + k3 * b_tableau[5][3] + k4 * b_tableau[5][4] + k5 * b_tableau[5][5], f)

    x4 = y + k1 * b_tableau[6][1] + k2 * b_tableau[6][2] + k3 * b_tableau[6][3] + k4 * b_tableau[6][4] + k5 * b_tableau[6][5]
    x5 = y + k1 * b_tableau[7][1] + k2 * b_tableau[7][2] + k3 * b_tableau[7][3] + k4 * b_tableau[7][4] + k5 * b_tableau[7][5] + k6 * b_tableau[7][6]
    return x4, x5

# Returns the state of the system relative to the solution for:
# some 4th order value "b4"
# some 5th order value "b5"
# where "0" = acceptable truncation error, "1" = truncation error too high, and "2" = truncation error too low
def check_step(b4, b5):
    if np.abs(b4 - b5) < truncation_error_upper:
        return True
    return False

# Returns the adjusted step size due to the amount of error in the previous step for:
# the previous step size "step_size"
# the 4th order value "b4" (used to calculate error)
# the 5th order value "b5" (used to calculate error)
def adjust_step(step_size, b4, b5):
    step_size = step_size * rk45_step_constant * (truncation_error_upper / np.abs(b4 - b5))**(1/4)
    return step_size

i = 0
x = 0 # initial x value
y = 0.5 # initial y value
step_size = 0.05 # initial step size
lambda_func = lambda x, y: y - x**2 + 1
while i < 20:
    errors = rk_45_step(x, y, step_size, lambda_func)
    val = check_step(errors[0], errors[1])
    if val:
        x += step_size
        y = errors[1]
        print(f"Step {i + 1}: {x}, {errors[1]} \nLocal Error: {np.abs(errors[0] - errors[1])}")
        step_size = adjust_step(step_size, errors[0], errors[1])
        i += 1
    else:
        step_size = adjust_step(step_size, errors[0], errors[1])
