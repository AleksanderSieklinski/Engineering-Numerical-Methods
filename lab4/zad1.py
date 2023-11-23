import numpy as np
import matplotlib.pyplot as plt
import concurrent.futures
import numba as nb

@nb.jit(nopython=True)
def calculate_S(V, rho,S=0):
    for i in range(0, nx):
        for j in range(0, ny):
            S += (delta**2) * (0.5*(((V[i+1,j]-V[i,j])/delta)**2) + 0.5*(((V[i,j+1]-V[i,j])/delta)**2) - (rho[i,j]*V[i,j]))
    return S
@nb.jit(nopython=True)
def global_relaxation(V, rho, omega):
    V_new = np.zeros_like(V)
    V_new[:, 0] = V1
    V_new[:, ny] = V2
    for i in range(1, nx):
        for j in range(1, ny):
            V_new[i, j] = 0.25 * (V[i+1, j] + V[i-1, j] + V[i, j+1] + V[i, j-1] + (delta**2/epsilon) * rho[i, j])
    for j in range (1, ny):
        V_new[0, j] = V_new[1, j]
        V_new[nx, j] = V_new[nx-1, j]
    V = (1 - omega) * V + omega * V_new
    return V
@nb.jit(nopython=True)
def local_relaxation(V, rho, omega):
    for i in range(1, nx):
        for j in range(1, ny):
            V[i, j] = (1 - omega) * V[i, j] + (omega/4) * (V[i+1, j] + V[i-1, j] + V[i, j+1] + V[i, j-1] + (delta**2/epsilon) * rho[i, j])
    for j in range (1, ny):
        V[0, j] = V[1, j]
        V[nx, j] = V[nx-1, j]
    return V
@nb.jit(nopython=True)
def rho(x, y):
    rho1 = np.exp(-((x - 0.35 * xmax)**2 / sigma_x**2 + (y - 0.5 * ymax)**2 / sigma_y**2))
    rho2 =-np.exp(-((x - 0.65 * xmax)**2 / sigma_x**2 + (y - 0.5 * ymax)**2 / sigma_y**2))
    return rho1 + rho2
def laplacian(V):
    grad_x, grad_y = np.gradient(V)
    grad_xx, _ = np.gradient(grad_x)
    _, grad_yy = np.gradient(grad_y)
    return (grad_xx + grad_yy)/(delta**2)
def solve_error(V, rho):
    delta = laplacian(V) + rho/epsilon
    # make all negative values positive
    delta[delta < 0] *= -1
    return delta
@nb.jit(nopython=True)
def compute_for_omega(omega, relaxation_type):
    V = np.zeros((nx+1, ny+1))
    V[:, 0] = V1
    V[:, ny] = V2
    S_list = []
    S = calculate_S(V, rho_array)
    it = 0
    if relaxation_type == 'global':
        while True:
            V = global_relaxation(V, rho_array, omega)
            it += 1
            S_new = calculate_S(V, rho_array)
            S_list.append(S_new)
            check = abs((S_new - S)/S)
            if check < TOL:
                break
            S = S_new
    elif relaxation_type == 'local':
        while True:
            V = local_relaxation(V, rho_array, omega)
            it += 1
            S_new = calculate_S(V, rho_array)
            S_list.append(S_new)
            check = abs((S_new - S)/S)
            if check < TOL:
                break
            S = S_new
    else:
        raise ValueError("Invalid relaxation_type. Expected 'global' or 'local'.")
    return S_list, it, V

epsilon = 1
delta = 0.1
nx = 150
ny = 100
V1 = 10
V2 = 0
xmax = delta * nx
ymax = delta * ny
sigma_x = 0.1 * xmax
sigma_y = 0.1 * ymax
omega_G = [0.6, 1.0]
omega_L = [1.0, 1.4, 1.8, 1.9]
TOL = 1e-8
delta_dict = {}
V_dict = {}

rho_array = np.zeros((nx+1, ny+1))
for i in range(nx+1):
    for j in range(ny+1):
        rho_array[i, j] = rho(i*delta, j*delta)

fig, (ax1, ax2) = plt.subplots(1, 2)

with concurrent.futures.ThreadPoolExecutor() as executor:
    futures = []

    for omega in omega_G:
        futures.append((omega,'global', executor.submit(compute_for_omega, omega, 'global')))
    for omega in omega_L:
        futures.append((omega,'local', executor.submit(compute_for_omega, omega, 'local')))

    for omega_val,relaxation_type, future in futures:
        S_list, it, V = future.result()
        ax = ax1 if relaxation_type == 'global' else ax2
        ax.plot(S_list, label=f'omega={omega_val}, {len(S_list)} it')
        print(f'type={relaxation_type}, omega={omega_val}, it={it}')
        if relaxation_type == 'global' and (omega_val == 0.6 or omega_val == 1.0):
            V_dict[omega_val] = V
            delta_dict[omega_val] = solve_error(V, rho_array)

ax1.set_title('Relaksacja globalna')
ax1.set_xlabel('Nr iteracji')
ax1.set_ylabel('S')
ax1.set_xticks([1,10, 100, 1000, 10000, 100000])
ax1.set_xscale('log')
ax1.legend()
ax2.set_title('Relaksacja lokalna')
ax2.set_xlabel('Nr iteracji')
ax2.set_ylabel('S')
ax2.set_xticks([1,10, 100, 1000, 10000, 100000])
ax2.set_xscale('log')
ax2.legend()
plt.savefig("z1_relaxation.png",bbox_inches='tight',transparent=True)
plt.show()
plt.clf()

plt.imshow(delta_dict[0.6], cmap='jet', extent=[0, xmax, 0, ymax])
plt.colorbar()
plt.title('Relaksacja globalna w=0.6')
plt.xlabel('x')
plt.ylabel('y')
plt.savefig("z1_error_06.png",bbox_inches='tight',transparent=True)
plt.show()
plt.clf()

plt.imshow(delta_dict[1.0], cmap='jet', extent=[0, xmax, 0, ymax])
plt.colorbar()
plt.xlabel('x')
plt.ylabel('y')
plt.title('Relaksacja globalna w=1.0')
plt.savefig("z1_error_10.png",bbox_inches='tight',transparent=True)
plt.show()
plt.clf()

plt.imshow(V_dict[0.6], cmap='jet', extent=[0, xmax, 0, ymax])
plt.colorbar()
plt.xlabel('x')
plt.ylabel('y')
plt.title('Relaksacja globalna w=0.6')
plt.savefig("z1_V_06.png",bbox_inches='tight',transparent=True)
plt.show()
plt.clf()

plt.imshow(V_dict[1.0], cmap='jet', extent=[0, xmax, 0, ymax])
plt.colorbar()
plt.xlabel('x')
plt.ylabel('y')
plt.title('Relaksacja globalna w=1.0')
plt.savefig("z1_V_10.png",bbox_inches='tight',transparent=True)
plt.show()
plt.clf()