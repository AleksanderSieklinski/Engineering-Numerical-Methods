import matplotlib.pyplot as plt
import numpy as np

def plot_temperature(it, filename,type):
    # Load data from file
    data = np.loadtxt(filename)
    # Reshape data into a 2D array
    x = data[:, 0]
    y = data[:, 1]
    T = data[:, 2]
    nx = np.unique(x).size
    ny = np.unique(y).size
    X = x.reshape(nx, ny)
    Y = y.reshape(nx, ny)
    T = T.reshape(nx, ny)

    # Create plot
    plt.figure(figsize=(6, 5))
    plt.contourf(X, Y, T, 50)
    plt.colorbar()
    plt.title(f'it = {it}')
    plt.xlabel('x')
    plt.ylabel('y')
    if it < 1000:
        if type == "t":
            plt.savefig(f'res/temperature_0{it}.png')
        else:
            plt.savefig(f'res/laplacian_0{it}.png')
    else:
        if type == "t":
            plt.savefig(f'res/temperature_{it}.png')
        else:
            plt.savefig(f'res/laplacian_{it}.png')
    plt.close()

def main():
    # Plot temperature maps at it = 100, 200, 500, 1000, 2000
    for it in [100, 200, 500, 1000, 2000]:
        if it < 1000:
            plot_temperature(it, f'out/temperature_0{it}.dat',"t")
            plot_temperature(it, f'out/laplacian_0{it}.dat',"l")
        else:
            plot_temperature(it, f'out/temperature_{it}.dat',"t")
            plot_temperature(it, f'out/laplacian_{it}.dat',"l")

if __name__ == '__main__':
    main()