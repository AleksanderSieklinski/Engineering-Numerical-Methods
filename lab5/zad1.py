import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import colormaps
import time

# parametry poczatkowe
delta = 0.2
nx = 128
ny = 128
xmax = delta * nx
ymax = delta * ny
TOL = 10**(-8)
kArray = [16,8,4,2,1]

def relaxate(V,k):
	if k == 1:
		for i in range(k,nx-k,k):
			for j in range(k,ny-k,k):
				V[i][j] = 0.25*(V[i-k][j]+V[i+k][j]+V[i][j-k]+V[i][j+k])
	else:
		for i in range(k,nx,k):
			for j in range(k,ny,k):
				try:
					V[i][j] = 0.25*(V[i-k][j]+V[i+k][j]+V[i][j-k]+V[i][j+k])
				except IndexError:
					V[i][j] = 0.25*(V[i-k][j]+V[i+k-1][j]+V[i][j-k]+V[i][j+k-1])

def S(V,k):
    sum = 0
    for i in range(0,nx-k,k):
        for j in range(0,ny-k,k):
            sum += (((k*delta)**2)/2)*( 
            (((V[i+k][j]-V[i][j])/(2*k*delta)) + ((V[i+k][j+k]-V[i][j+k])/(2*k*delta)) )**2 + 
            (((V[i][j+k]-V[i][j])/(2*k*delta)) + ((V[i+k][j+k]-V[i+k][j])/(2*k*delta)) )**2 )
    return sum

def dens_mesh(V,k):
    k2 = k//2
    for i in range(0,nx-k,k):
        for j in range(k,ny-k,k):
            V[i+k2][j+k2] = 0.25*(V[i][j]+V[i+k][j]+V[i][j+k]+V[i+k][j+k])
            V[i+k][j+k2] = 0.5*(V[i+k][j]+V[i+k][j+k])
            V[i+k2][j+k] = 0.5*(V[i][j+k]+V[i+k][j+k])
            V[i+k2][j] = 0.5*(V[i][j]+V[i+k][j])
            V[i][j+k2] = 0.5*(V[i][j]+V[i][j+k])

def main():
    V = np.zeros((nx, ny))
    S_values = []
    sArray = {}
    iterArray = {}
    # warunki brzegowe
    for i in range(0, nx):
        V[0][i] = math.sin(math.pi*(i*delta/ymax)) #Vb1
        V[ny-1][i] = math.sin(math.pi*(i*delta/ymax)) #Vb3
    for j in range(0, ny):
        V[j][ny-1] = -math.sin(2*math.pi*(j*delta/xmax)) #Vb2
        V[j][0] = math.sin(2*math.pi*(j*delta/xmax)) #Vb4
    # zad2 mapy potencjału dla wszystkich k
    for k in kArray:
        print("k = ", k)
        sArr = []
        iterArr = []
        iter = 0
        while True:
            iterArr.append(iter)
            iter += 1
            V_old = V.copy()
            relaxate(V, k)
            S_value = S(V, k)
            print(S_value)
            S_values.append(S_value)
            sArr.append(S_value)
            if np.abs(S_value - S(V_old, k)) < TOL:
                break
        sArray[k] = sArr
        iterArray[k] = iterArr
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('k='+str(k))
        plt.imshow(np.rot90(V[::k, ::k]), cmap='turbo')
        plt.colorbar()
        plt.savefig('z1_k='+str(k)+'.png',bbox_inches='tight',transparent=True)
        plt.clf() 
        if k != 1:
            dens_mesh(V, k)
    # zad2 zmiany S(it) dla wszystkich k
    for k in kArray:
        plt.plot(iterArray[k],sArray[k], label='k='+str(k))
    plt.title("S(it)")
    plt.xlabel('it')
    plt.ylabel('S')
    plt.legend()
    plt.savefig('z1_S(it).png',bbox_inches='tight',transparent=True)
    plt.clf()

if __name__ == "__main__":
    main()