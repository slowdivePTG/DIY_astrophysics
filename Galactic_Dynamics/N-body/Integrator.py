import numpy as np

# Integrators solving 2nd order ODE: d^2x(t)/dt^2 = f(x(t))
# Initial condition: X = [x(t), dx(t)/dt]^T


def Runge_Kutta(X, t, f, dT):
    '''
    4th order Runge-Kutta integrator
    '''
    def F(X, t, f):
        dXdt1 = X[1]
        dXdt2 = f(X[0], t)
        return np.array([dXdt1, dXdt2])

    k1 = F(X, t, f)
    k2 = F(X + k1 * dT / 2, t + dT / 2, f)
    k3 = F(X + k2 * dT / 2, t + dT / 2, f)
    k4 = F(X + k3 * dT, t + dT, f)
    t += dT
    return X + (k1 / 6 + k2 / 3 + k3 / 3 + k4 / 6) * dT, t + dT


def LeapFrog(X, t, f, dT):
    '''
    Leapfrog integrator
    '''
    vtemp = X[1] + f(X[0], t) * dT / 2
    X[0] = vtemp * dT + X[0]
    X[1] = vtemp + f(X[0], t) * dT / 2
    t += dT
    return X, t
