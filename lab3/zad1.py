import numpy as np
import matplotlib.pyplot as plt

x0 = 0.01
v0 = 0
dt0 = 1
S = 0.75
tol1 = 1e-2
tol2 = 1e-5
tmax = 40
alpha = 5

def f(t, x, v, alpha):
    return v

def g(t, x, v, alpha):
    return alpha*(1-x**2)*v - x

def trapezoidal(xn, vn, dt, alpha):
    def F(xnp1, vnp1):
        return xnp1 - xn - 0.5*dt*(f(0, xn, vn, alpha) + f(0, xnp1, vnp1, alpha))
    
    def G(xnp1, vnp1):
        return vnp1 - vn - 0.5*dt*(g(0, xn, vn, alpha) + g(0, xnp1, vnp1, alpha))
    
    xnp1 = xn
    vnp1 = vn
    a11 = 1
    a12 = -0.5*dt
    while True:
        F_val = F(xnp1, vnp1)
        G_val = G(xnp1, vnp1)

        a21 = -0.5*dt*(-2*alpha*xnp1*vnp1 - 1)
        a22 = 1 - 0.5*dt*alpha*(1 - xnp1**2)

        dx = (-F_val*a22 + G_val*a12)/(a11*a22 - a12*a21)
        dv = (F_val*a21 - G_val*a11)/(a11*a22 - a12*a21)

        xnp1 += dx
        vnp1 += dv

        if max(abs(dx), abs(dv)) < 1e-10:
            break
    
    return xnp1, vnp1


def RK2(xn, vn, dt, alpha):
    k1x = f(0, xn, vn, alpha)
    k1v = g(0, xn, vn, alpha)
    k2x = f(0, xn + dt*k1x, vn + dt*k1v, alpha)
    k2v = g(0, xn + dt*k1x, vn + dt*k1v, alpha)
    xn1 = xn + 0.5*dt*(k1x + k2x)
    vn1 = vn + 0.5*dt*(k1v + k2v)
    return xn1, vn1

def control_step_size(dt, S, tol, p, Ex, Ev):
    return ((S*tol)/max(abs(Ex), abs(Ev)))**(1/(p+1))*dt

def solve_equation(method, xn, vn, dt, alpha, tol, tmax, S=0.75, p=2):
    t = 0
    x_list = [xn]
    v_list = [vn]
    dt_list = [dt]
    t_list = [t]
    while t < tmax:
        xn1, vn1 = method(xn, vn, dt, alpha)
        xn2, vn2 = method(xn1, vn1, dt, alpha)
        xn21, vn21 = method(xn, vn, 2*dt, alpha)
        Ex = (xn2 - xn21)/(2**p - 1)
        Ev = (vn2 - vn21)/(2**p - 1)
        if max(abs(Ex), abs(Ev)) < tol:
            t += 2*dt
            xn = xn2
            vn = vn2
            x_list.append(xn)
            v_list.append(vn)
            dt_list.append(dt)
            t_list.append(t)
        dt = control_step_size(dt,S,tol,p,Ex,Ev)
        print(t, dt,tmax)
    return x_list, v_list, dt_list,t_list

x_list1, v_list1, dt_list1,t_list1 = solve_equation(RK2, x0, v0, dt0, alpha, tol1, tmax)
x_list2, v_list2, dt_list2,t_list2 = solve_equation(RK2, x0, v0, dt0, alpha, tol2, tmax)

plt.figure(figsize=(10, 8))
plt.subplot(2, 2, 1)
plt.plot(t_list1, x_list1)
plt.plot(t_list2, x_list2)
plt.xlim([0, tmax])
plt.title('Metoda RK2')
plt.xlabel('t')
plt.ylabel('x(t)')
plt.legend(['TOL=1e-2', 'TOL=1e-5'])
plt.subplot(2, 2, 2)
plt.plot(t_list1, v_list1)
plt.plot(t_list2, v_list2)
plt.xlim([0, tmax])
plt.title('Metoda RK2')
plt.xlabel('t')
plt.ylabel('v(t)')
plt.legend(['TOL=1e-2', 'TOL=1e-5'])
plt.subplot(2, 2, 3)
plt.plot(t_list1, dt_list1, 'o-', markersize=3)  
plt.plot(t_list2, dt_list2, 'o-', markersize=3)
plt.xlim([0, tmax])
plt.title('Metoda RK2')
plt.xlabel('t')
plt.ylabel('dt(t)')
plt.legend(['TOL=1e-2', 'TOL=1e-5'])
plt.subplot(2, 2, 4)
plt.plot(x_list1, v_list1)
plt.plot(x_list2, v_list2)
plt.title('Metoda RK2')
plt.xlabel('x')
plt.ylabel('v')
plt.legend(['TOL=1e-2', 'TOL=1e-5'])
plt.show()

x_list1, v_list1, dt_list1,t_list1 = solve_equation(trapezoidal, x0, v0, dt0, alpha, tol1, tmax)
x_list2, v_list2, dt_list2,t_list2 = solve_equation(trapezoidal, x0, v0, dt0, alpha, tol2, tmax)

plt.figure(figsize=(10, 8))
plt.subplot(2, 2, 1)
plt.plot(t_list1, x_list1)  
plt.plot(t_list2, x_list2)
plt.xlim([0, tmax])
plt.title('Metoda trapezow')
plt.xlabel('t')
plt.ylabel('x(t)')
plt.legend(['TOL=1e-2', 'TOL=1e-5'])
plt.subplot(2, 2, 2)
plt.plot(t_list1, v_list1)  
plt.plot(t_list2, v_list2)
plt.xlim([0, tmax])
plt.title('Metoda trapezow')
plt.xlabel('t')
plt.ylabel('v(t)')
plt.legend(['TOL=1e-2', 'TOL=1e-5'])
plt.subplot(2, 2, 3)
plt.plot(t_list1, dt_list1, 'o-', markersize=3)  
plt.plot(t_list2, dt_list2, 'o-', markersize=3)
plt.xlim([0, tmax])
plt.title('Metoda trapezow')
plt.xlabel('t')
plt.ylabel('dt(t)')
plt.legend(['TOL=1e-2', 'TOL=1e-5'])
plt.subplot(2, 2, 4)
plt.plot(x_list1, v_list1)
plt.plot(x_list2, v_list2)
plt.title('Metoda trapezow')
plt.xlabel('x')
plt.ylabel('v')
plt.legend(['TOL=1e-2', 'TOL=1e-5'])
plt.show()