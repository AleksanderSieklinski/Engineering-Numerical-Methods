import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def euler(y0, lambda_, t0, tn, dt):
    t = np.arange(t0, tn, dt)
    y = np.zeros(len(t))
    y[0] = y0
    for i in range(1, len(t)):
        y[i] = y[i-1] + dt*lambda_*y[i-1]
    return t, y

def rk2(y0, lambda_, t0, tn, dt):
    t = np.arange(t0, tn, dt)
    y = np.zeros(len(t))
    y[0] = y0
    for i in range(1, len(t)):
        k1 = lambda_*y[i-1]
        k2 = lambda_*(y[i-1] + dt*k1)
        y[i] = y[i-1] + dt/2*(k1 + k2)
    return t, y

def rk4(y0, lambda_, t0, tn, dt):
    t = np.arange(t0, tn, dt)
    y = np.zeros(len(t))
    y[0] = y0
    for i in range(1, len(t)):
        k1 = lambda_*y[i-1]
        k2 = lambda_*(y[i-1] + dt/2*k1)
        k3 = lambda_*(y[i-1] + dt/2*k2)
        k4 = lambda_*(y[i-1] + dt*k3)
        y[i] = y[i-1] + dt/6*(k1 + 2*k2 + 2*k3 + k4)
    return t, y

def euler_global_error(y0, lambda_, t0, tn, dt):
    t, y = euler(y0, lambda_, t0, tn, dt)
    y_analytical = np.exp(lambda_*t)
    return y - y_analytical

def rk2_global_error(y0, lambda_, t0, tn, dt):
    t, y = rk2(y0, lambda_, t0, tn, dt)
    y_analytical = np.exp(lambda_*t)
    return y - y_analytical

def rk4_global_error(y0, lambda_, t0, tn, dt):
    t, y = rk4(y0, lambda_, t0, tn, dt)
    y_analytical = np.exp(lambda_*t)
    return y - y_analytical

def f2(Q, I, t, R, L, C, V0, omega):
    return I

def g2(Q, I, t, R, L, C, V0, omega):
    return V0*np.sin(omega*t)/L - R/L*I - 1/(L*C)*Q

def rk4_2(Q0, I0, t0, tn, dt, R, L, C, V0, omega):
    t = np.arange(t0, tn, dt)
    Q = np.zeros(len(t))
    I = np.zeros(len(t))
    Q[0] = Q0
    I[0] = I0
    for i in range(1, len(t)):
        kQ1 = f2(Q[i-1], I[i-1], t[i-1], R, L, C, V0, omega)
        kI1 = g2(Q[i-1], I[i-1], t[i-1], R, L, C, V0, omega)
        kQ2 = f2(Q[i-1] + dt/2*kQ1, I[i-1] + dt/2*kI1, t[i-1] + dt/2, R, L, C, V0, omega)
        kI2 = g2(Q[i-1] + dt/2*kQ1, I[i-1] + dt/2*kI1, t[i-1] + dt/2, R, L, C, V0, omega)
        kQ3 = f2(Q[i-1] + dt/2*kQ2, I[i-1] + dt/2*kI2, t[i-1] + dt/2, R, L, C, V0, omega)
        kI3 = g2(Q[i-1] + dt/2*kQ2, I[i-1] + dt/2*kI2, t[i-1] + dt/2, R, L, C, V0, omega)
        kQ4 = f2(Q[i-1] + dt*kQ3, I[i-1] + dt*kI3, t[i-1] + dt, R, L, C, V0, omega)
        kI4 = g2(Q[i-1] + dt*kQ3, I[i-1] + dt*kI3, t[i-1] + dt, R, L, C, V0, omega)
        Q[i] = Q[i-1] + dt/6*(kQ1 + 2*kQ2 + 2*kQ3 + kQ4)
        I[i] = I[i-1] + dt/6*(kI1 + 2*kI2 + 2*kI3 + kI4)
    return t, Q, I

def main():
    y0 = 1
    lambda_ = -1
    t0 = 0
    tn = 5
    dt = [0.01, 0.1, 1.0]


    # zad1 rozw
    t, y = euler(y0, lambda_, t0, tn, dt[0])
    plt.plot(t, y, label=f"dt={dt[0]}", linestyle='None', marker='o')
    t, y = euler(y0, lambda_, t0, tn, dt[1])
    plt.plot(t, y, label=f"dt={dt[1]}", linestyle = 'dashed')
    t, y = euler(y0, lambda_, t0, tn, dt[2])
    plt.plot(t, y, label=f"dt={dt[2]}", linestyle = 'none', marker='o')
    t = np.arange(t0, tn, 0.01)
    y = np.exp(lambda_*t)
    plt.plot(t, y, label="analityczne")
    plt.title("Euler rozwiazanie")
    plt.xlabel("t")
    plt.ylabel("y(t)")
    plt.legend()
    plt.savefig("z1_metoda_eulera_rozw.png",bbox_inches='tight',transparent=True)
    plt.show()
    plt.clf()
    # zad1 blad
    abs = euler_global_error(y0, lambda_, t0, tn, dt[0])
    t = np.arange(t0, tn, dt[0])
    plt.plot(t, abs, label=f"dt={dt[0]}")
    abs = euler_global_error(y0, lambda_, t0, tn, dt[1])
    t = np.arange(t0, tn, dt[1])
    plt.plot(t, abs, label=f"dt={dt[1]}")
    abs = euler_global_error(y0, lambda_, t0, tn, dt[2])
    t = np.arange(t0, tn, dt[2])
    plt.plot(t, abs, label=f"dt={dt[2]}")
    plt.xlabel("t")
    plt.ylabel("ynum(t) - ydok(t)")
    plt.title("Euler błąd globalny")
    plt.legend()
    plt.savefig("z1_metoda_eulera_blad.png",bbox_inches='tight',transparent=True)
    plt.show()
    plt.clf()
    # zad2 rozw
    t, y = rk2(y0, lambda_, t0, tn, dt[0])
    plt.plot(t, y, label=f"dt={dt[0]}", linestyle='None', marker='o')
    t, y = rk2(y0, lambda_, t0, tn, dt[1])
    plt.plot(t, y, label=f"dt={dt[1]}", linestyle = 'dashed')
    t, y = rk2(y0, lambda_, t0, tn, dt[2])
    plt.plot(t, y, label=f"dt={dt[2]}", linestyle = 'none', marker='o')
    t = np.arange(t0, tn, 0.01)
    y = np.exp(lambda_*t)
    plt.plot(t, y, label="analityczne")
    plt.title("RK2 rozwiazanie")
    plt.xlabel("t")
    plt.ylabel("y(t)")
    plt.legend()
    plt.savefig("z2_rk2_rozw.png",bbox_inches='tight',transparent=True)
    plt.show()

    plt.clf()
    # zad2 blad
    abs = rk2_global_error(y0, lambda_, t0, tn, dt[0])
    t = np.arange(t0, tn, dt[0])
    plt.plot(t, abs, label=f"dt={dt[0]}")
    abs = rk2_global_error(y0, lambda_, t0, tn, dt[1])
    t = np.arange(t0, tn, dt[1])
    plt.plot(t, abs, label=f"dt={dt[1]}")
    abs = rk2_global_error(y0, lambda_, t0, tn, dt[2])
    t = np.arange(t0, tn, dt[2])
    plt.plot(t, abs, label=f"dt={dt[2]}")
    plt.xlabel("t")
    plt.ylabel("ynum(t) - ydok(t)")
    plt.title("RK2 błąd globalny")
    plt.legend()
    plt.savefig("z2_rk2_blad.png",bbox_inches='tight',transparent=True)
    plt.show()
    plt.clf()
    # zad3 rozw
    t, y = rk4(y0, lambda_, t0, tn, dt[0])
    plt.plot(t, y, label=f"dt={dt[0]}", linestyle='None', marker='o')
    t, y = rk4(y0, lambda_, t0, tn, dt[1])
    plt.plot(t, y, label=f"dt={dt[1]}", linestyle = 'dashed')
    t, y = rk4(y0, lambda_, t0, tn, dt[2])
    plt.plot(t, y, label=f"dt={dt[2]}", linestyle = 'none', marker='o')
    t = np.arange(t0, tn, 0.01)
    y = np.exp(lambda_*t)
    plt.plot(t, y, label="analityczne")
    plt.title("RK4 rozwiazanie")
    plt.xlabel("t")
    plt.ylabel("y(t)")
    plt.legend()
    plt.savefig("z3_rk4_rozw.png",bbox_inches='tight',transparent=True)
    plt.show()
    plt.clf()
    # zad3 blad
    abs = rk4_global_error(y0, lambda_, t0, tn, dt[0])
    t = np.arange(t0, tn, dt[0])
    plt.plot(t, abs, label=f"RK4 dt={dt[0]}")
    abs = rk4_global_error(y0, lambda_, t0, tn, dt[1])
    t = np.arange(t0, tn, dt[1])
    plt.plot(t, abs, label=f"RK4 dt={dt[1]}")
    abs = rk4_global_error(y0, lambda_, t0, tn, dt[2])
    t = np.arange(t0, tn, dt[2])
    plt.plot(t, abs, label=f"RK4 dt={dt[2]}")
    plt.xlabel("t")
    plt.ylabel("ynum(t) - ydok(t)")
    plt.title("RK4 błąd globalny")
    plt.legend()
    plt.savefig("z3_rk4_blad.png",bbox_inches='tight',transparent=True)
    plt.show()
    plt.clf()


    Q0 = 0
    I0 = 0
    t0 = 0
    tn = 4*2*np.pi/np.sqrt(1/(0.1*0.001))
    dt = 10**(-4)
    R = 100
    L = 0.1
    C = 0.001
    V0 = 10
    # zad4 rozw Q(t)
    omega = [0.5*np.sqrt(1/(0.1*0.001)), 0.8*np.sqrt(1/(0.1*0.001)), 1.0*np.sqrt(1/(0.1*0.001)), 1.2*np.sqrt(1/(0.1*0.001))]
    t, Q, I = rk4_2(Q0, I0, t0, tn, dt, R, L, C, V0, omega[0])
    plt.plot(t, Q, label=f"0.5*w0")
    t, Q, I = rk4_2(Q0, I0, t0, tn, dt, R, L, C, V0, omega[1])
    plt.plot(t, Q, label=f"0.8*w0")
    t, Q, I = rk4_2(Q0, I0, t0, tn, dt, R, L, C, V0, omega[2])
    plt.plot(t, Q, label=f"1.0*w0")
    t, Q, I = rk4_2(Q0, I0, t0, tn, dt, R, L, C, V0, omega[3])
    plt.plot(t, Q, label=f"1.2*w0")
    plt.xlabel("t")
    plt.ylabel("Q(t)")
    plt.title("Metoda RK4 Q(t)")
    plt.legend()
    plt.savefig("z4_rk4_Q(t).png",bbox_inches='tight',transparent=True)
    plt.show()
    plt.clf()
    # zad4 rozw I(t)
    t, Q, I = rk4_2(Q0, I0, t0, tn, dt, R, L, C, V0, omega[0])
    plt.plot(t, I,label=f"0.5*w0")
    t, Q, I = rk4_2(Q0, I0, t0, tn, dt, R, L, C, V0, omega[1])
    plt.plot(t, I,label=f"0.8*w0")
    t, Q, I = rk4_2(Q0, I0, t0, tn, dt, R, L, C, V0, omega[2])
    plt.plot(t, I,label=f"1.0*w0")
    t, Q, I = rk4_2(Q0, I0, t0, tn, dt, R, L, C, V0, omega[3])
    plt.plot(t, I,label=f"1.2*w0")
    plt.xlabel("t")
    plt.ylabel("I(t)")
    plt.title("Metoda RK4 I(t)")
    plt.legend()
    plt.savefig("z4_rk4_I(t).png",bbox_inches='tight',transparent=True)
    plt.show()
    plt.clf()


if __name__ == "__main__":
    main()