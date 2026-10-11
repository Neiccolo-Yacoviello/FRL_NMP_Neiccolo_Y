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

# def check_step():


print(rk_45_step(0, 0.5, 0.2, lambda x, y: y - x**2 + 1))