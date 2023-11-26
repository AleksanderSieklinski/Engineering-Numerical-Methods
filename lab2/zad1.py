import numpy as np
import matplotlib.pyplot as plt

beta = 0.001
N = 500
gamma = 0.1
tmax = 100
dt = 0.1
u0 = 1
TOL = 1e-6
max_iter = 20

def f(t, u):
    alpha = beta * N - gamma
    return (alpha * u - beta * u**2)

def trapezoidal_picard():
    t = np.arange(0, tmax + dt, dt)
    u = np.zeros(len(t))
    u[0] = u0
    z = N - u
    for i in range(1, len(t)):
        un = u[i-1]
        for j in range(max_iter):
            u[i] = u[i-1] + dt/2 * (f(t[i], un) + f(t[i-1], u[i-1]))
            if abs(u[i] - un) < TOL:
                break
            un = u[i]
        z[i] = N - u[i]
    return t, u, z

def trapezoidal_newton():
    t = np.arange(0, tmax + dt, dt)
    u = np.zeros(len(t))
    u[0] = u0
    z = N - u
    for i in range(1, len(t)):
        un = u[i-1]
        for j in range(max_iter):
            F = lambda x: x - u[i-1] - dt/2 * (f(t[i-1], u[i-1]) + f(t[i], x))
            dF = lambda x: 1 - dt/2 * beta * (N - 2 * x)
            u[i] = un - F(un) / dF(un)
            if abs(u[i] - un) < TOL:
                break
            un = u[i]
        z[i] = N - u[i]
    return t, u, z

def RK2(tmax=tmax, dt=dt, u0=u0, TOL=TOL, max_iter=max_iter):
    t = np.arange(0, tmax + dt, dt)
    u = np.zeros(len(t))
    u[0] = u0
    z = N - u
    for i in range(1, len(t)):
        un = u[i-1]
        for j in range(max_iter):
            F1 = lambda x1, x2: x1 - un - dt * (a11 * f(t[i-1] + c1*dt, x1) + a12 * f(t[i-1] + c2*dt, x2))
            F2 = lambda x1, x2: x2 - un - dt* (a21 * f(t[i-1] + c1*dt, x1) + a22 * f(t[i-1] + c2*dt, x2))
            dF11 = lambda x1, x2: 1 - dt * a11 * beta * (N - 2 * x1)
            dF12 = lambda x1, x2: -dt * a12 * beta * (N - 2 * x2)
            dF21 = lambda x1, x2: -dt * a21 * beta * (N - 2 * x1)
            dF22 = lambda x1, x2: 1 - dt * a22 * beta * (N - 2 * x2)
            m11 = dF11(un, un)
            m12 = dF12(un, un)
            m21 = dF21(un, un)
            m22 = dF22(un, un)
            det = m11 * m22 - m12 * m21
            if det == 0:
                break
            inv_m11 = m22 / det
            inv_m12 = -m12 / det
            inv_m21 = -m21 / det
            inv_m22 = m11 / det
            U1 = un
            U2 = un
            for k in range(max_iter):
                dU1 = inv_m11 * F1(U1, U2) + inv_m12 * F2(U1, U2)
                dU2 = inv_m21 * F1(U1, U2) + inv_m22 * F2(U1, U2)
                U1 -= dU1
                U2 -= dU2
                if abs(dU1) < TOL and abs(dU2) < TOL:
                    break
            u[i] = un + dt * (b1 * f(t[i-1] + c1*dt, U1) + b2 * f(t[i-1] + c2*dt, U2))
            z[i] = N - u[i]
    return t, u, z

# zad1 Picard
t1, u1, z1 = trapezoidal_picard()
t2, u2, z2 = trapezoidal_newton()
plt.plot(t1, u1, label='u(t)')
plt.plot(t1, z1, label='v(t)')
plt.title('Metoda Picarda')
plt.xlabel('t')
plt.ylabel('u(t),v(t)')
plt.xticks(np.arange(0, tmax + dt, 50))
plt.yticks(np.arange(0, N + 1, 200))
plt.xlim([0, t1[-1]])
plt.ylim([0, N])
plt.legend()
plt.savefig("z1_Picard.png",bbox_inches='tight',transparent=True)
plt.show()
plt.clf()
# zad1 Newton
plt.plot(t2, u2, label='u(t)')
plt.plot(t2, z2, label='v(t)')
plt.title('Metoda Newtona')
plt.xlabel('t')
plt.ylabel('u(t),v(t)')
plt.xticks(np.arange(0, tmax + dt, 50))
plt.yticks(np.arange(0, N + 1, 200))
plt.xlim([0, t2[-1]])
plt.ylim([0, N])
plt.legend()
plt.savefig("z1_Newton.png",bbox_inches='tight',transparent=True)
plt.show()
plt.clf()
# zad2 RK2
c1 = 1/2 - np.sqrt(3)/6
c2 = 1/2 + np.sqrt(3)/6
a11 = a22 = 1/4
a12 = a21 = 1/4 - np.sqrt(3)/6
b1 = 1/2 + np.sqrt(3)/6
b2 = 1/2 - np.sqrt(3)/6
t, u, z = RK2()
plt.plot(t, u, label='u(t)')
plt.plot(t, z, label='v(t)')
plt.title('Niejawna RK2')
plt.xlabel('Time')
plt.ylabel('Population(u(t),v(t))')
plt.xticks(np.arange(0, tmax + dt, 50))
plt.yticks(np.arange(0, N + 1, 200))
plt.xlim([0, t[-1]])
plt.ylim([0, N])
plt.legend()
plt.savefig("z2_RK2.png",bbox_inches='tight',transparent=True)
plt.show()
plt.clf()