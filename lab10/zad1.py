import numpy as np
import matplotlib.pyplot as plt

def plot_energy(a, b, ax=None):
    E = np.loadtxt('out/energy_alfa_{0}_beta_{1}.txt'.format(a, b))
    t = np.linspace(0, 50, len(E))
    if ax is None:
        ax = plt.gca()
    ax.plot(t, E, label='a={0}, b={1}'.format(a, b))
    ax.set_xlabel('t')
    ax.set_ylabel('E')
    plt.xticks(np.arange(0, 51, 5))
    ax.set_title('E(t)')
    ax.legend()
# wykres transponowany, bo w pliku zapisane są wierszami a nie kolumnami
def plot_wave(a,b):
    u = np.loadtxt('out/output_alfa_{0}_beta_{1}.txt'.format(a, b)).T
    x = np.linspace(0, 16, u.shape[1])
    t = np.linspace(0, 50, u.shape[0])
    plt.figure(figsize=(10, 5))
    plt.imshow(u, aspect='auto', cmap='hot', extent=[t[0], t[-1], x[0], x[-1]])
    plt.colorbar()
    plt.xlabel('t')
    plt.ylabel('x')
    plt.xticks(np.arange(0, 51, 5))
    plt.title('B={0}, a={1}'.format(b, a))
    plt.savefig('res/wave_alfa{0}_beta{1}.png'.format(a,b))

def main():
    a = [0.0,0.0,0.0,1.0]
    b = [0.0,0.1,1.0,1.0]
    fig, ax = plt.subplots(figsize=(10, 5))
    for i in range(3):
        plot_energy(a[i], b[i], ax)
    ax.figure.savefig('res/energy.png')
    plt.clf()
    plot_energy(a[3], b[3])
    plt.savefig('res/energy_alfa_1.0_beta_1.0.png')
    for i in range(len(b)):
        plot_wave(a[i], b[i])

if __name__ == '__main__':
    main()