import numpy as np
import matplotlib.pyplot as plt
import math
from matplotlib import colormaps
import numba as nb

delta = 0.2
nx = 128
ny = 128
xmax = delta * nx
ymax = delta * ny
TOL = 1e-8
kArray = [16,8,4,2,1]

@nb.jit(nopython=True)
def relaxate(V,k):
    for i in range(k,nx,k):
        for j in range(k,ny,k):
            V[i][j] = 0.25*(V[i+k][j] + V[i-k][j] + V[i][j+k] + V[i][j-k])
@nb.jit(nopython=True)
def S(V,k):
    sum = 0
    for i in range(0,nx,k):
        for j in range(0,ny,k):
            sum += 0.5 * pow(k*delta,2)*(
            pow((V[i+k,j] - V[i,j] + V[i+k,j+k] - V[i,j+k])/(2*k*delta),2)+
            pow((V[i,j+k] - V[i,j] + V[i+k,j+k] - V[i+k,j])/(2*k*delta),2) 
            )
    return sum
@nb.jit(nopython=True)
def dens_mesh(V,k):
    k2 = k//2
    for i in range(0,nx,k):
        for j in range(0,ny,k):
            V[i+k2][j+k2] = 0.25*(V[i][j]+V[i+k][j]+V[i][j+k]+V[i+k][j+k])
            V[i+k][j+k2] = 0.5*(V[i+k][j]+V[i+k][j+k])
            V[i+k2][j+k] = 0.5*(V[i][j+k]+V[i+k][j+k])
            V[i+k2][j] = 0.5*(V[i][j]+V[i+k][j])
            V[i][j+k2] = 0.5*(V[i][j]+V[i][j+k])

def main():
    V = np.zeros((nx+1, ny+1))
    sArray = {}
    iterArray = {}
    # zad2 mapy potencjału dla wszystkich k
    iter = 0
    for k in kArray:
        # warunki brzegowe
        for i in range(0, nx):
            V[i][ny] = -math.sin(2*math.pi*(i*delta/xmax)) #Vb2
            V[i][0] = math.sin(2*math.pi*(i*delta/xmax)) #Vb4
        for j in range(0, ny):
            V[0][j] = math.sin(math.pi*(j*delta/ymax)) #Vb1
            V[nx][j] = math.sin(math.pi*(j*delta/ymax)) #Vb3
        print("k = ", k)
        sArr = []
        iterarr = []
        s_prev = 1
        while True:
            iter += 1
            iterarr.append(iter)
            relaxate(V, k)
            S_value = S(V, k)
            sArr.append(S_value)
            if np.abs((S_value - s_prev)/s_prev) < TOL:
                break
            s_prev = S_value
        sArray[k] = sArr
        iterArray[k] = iterarr
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('k='+str(k))
        plt.imshow(np.rot90(V[::k, ::k]), cmap='turbo',vmin=-1,vmax=1, extent=[0, xmax, 0, ymax])
        plt.xticks(np.arange(0, 26, 5))
        plt.yticks(np.arange(0, 26, 5))
        plt.colorbar(ticks=[-1, -0.5, 0, 0.5, 1])
        plt.savefig('z1_k='+str(k)+'.png',bbox_inches='tight',transparent=True)
        plt.clf() 
        if k != 1:
            dens_mesh(V, k)
    # zad2 zmiany S(it) dla wszystkich k
    for k in kArray:
        plt.plot(iterArray[k], sArray[k], label='k='+str(k))
    plt.title("S(it)")
    plt.xlabel('it')
    plt.ylabel('S')
    plt.legend()
    plt.xticks(np.arange(0, 700, 100))
    plt.yticks(np.arange(4.2, 5.65, 0.2))
    plt.savefig('z1_S(it).png',bbox_inches='tight',transparent=True)
    plt.clf()

if __name__ == "__main__":
    main()