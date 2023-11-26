import matplotlib.pyplot as plt
import numpy as np
# zad 5 podpunkt a
arr = []
with open("z5_nx_ny_50.txt") as file:
	for line in file.readlines():
			arr.append([float(i) for i in line.split(' ') if i != '\n'])
fig, ax = plt.subplots()
plt.xlabel('x')
plt.ylabel('y')
plt.title('nx=ny=50, e1=1,e2=1')
x_max = 50
y_max = 50
plt.xticks(np.arange(0, x_max/10 + 1, 5))
plt.yticks(np.arange(0, y_max/10 + 1, 5))
plt.imshow(arr, cmap='seismic',vmin=-10,vmax=10, extent=[0, x_max/10, 0, y_max/10])
plt.colorbar(ticks=[-10,0,10])
plt.savefig('z5_nx_ny_50_plot.png')
# zad 5 podpunkt b
arr = []
with open("z5_nx_ny_100.txt") as file:
	for line in file.readlines():
			arr.append([float(i) for i in line.split(' ') if i != '\n'])
fig, ax = plt.subplots()
plt.xlabel('x')
plt.ylabel('y')
plt.title('nx=ny=100')
x_max = 100
y_max = 100
plt.xticks(np.arange(0, x_max/10 + 1, 5))
plt.yticks(np.arange(0, y_max/10 + 1, 5))
plt.imshow(arr, cmap='seismic',vmin=-10,vmax=10, extent=[0, x_max/10, 0, y_max/10])
plt.colorbar(ticks=[-10,0,10])
plt.savefig('z5_nx_ny_100_plot.png')
# zad 5 podpunkt c
arr = []
with open("z5_nx_ny_200.txt") as file:
	for line in file.readlines():
			arr.append([float(i) for i in line.split(' ') if i != '\n'])
fig, ax = plt.subplots()
plt.xlabel('x')
plt.ylabel('y')
plt.title('nx=ny=200')
x_max = 200
y_max = 200
plt.xticks(np.arange(0, x_max/10 + 1, 5))
plt.yticks(np.arange(0, y_max/10 + 1, 5))
plt.imshow(arr, cmap='seismic',vmin=-10,vmax=10, extent=[0, x_max/10, 0, y_max/10])
plt.colorbar(ticks=[-10,0,10])
plt.savefig('z5_nx_ny_200_plot.png')
# zad 6 podpunkt a
arr = []
with open("z6_e1_1_e2_1.txt") as file:
	for line in file.readlines():
			arr.append([float(i) for i in line.split(' ') if i != '\n'])
fig, ax = plt.subplots()
plt.xlabel('x')
plt.ylabel('y')
plt.title('nx=ny=100, e1=1,e2=1')
x_max = 100
y_max = 100
plt.xticks(np.arange(0, x_max/10 + 1, 5))
plt.yticks(np.arange(0, y_max/10 + 1, 5))
plt.imshow(arr, cmap='seismic',vmin=-1,vmax=1, extent=[0, x_max/10, 0, y_max/10])
plt.colorbar(ticks=[-0.8,0,0.8])
plt.savefig('z6_e1_1_e2_1_plot.png')
# zad 6 podpunkt b
arr = []
with open("z6_e1_1_e2_2.txt") as file:
	for line in file.readlines():
			arr.append([float(i) for i in line.split(' ') if i != '\n'])
fig, ax = plt.subplots()
plt.xlabel('x')
plt.ylabel('y')
plt.title('e1=1,e2=2')
x_max = 100
y_max = 100
plt.xticks(np.arange(0, x_max/10 + 1, 5))
plt.yticks(np.arange(0, y_max/10 + 1, 5))
plt.imshow(arr, cmap='seismic',vmin=-1,vmax=1, extent=[0, x_max/10, 0, y_max/10])
plt.colorbar(ticks=[-0.8,0,0.8])
plt.savefig('z6_e1_1_e2_2_plot.png')
# zad 6 podpunkt c
arr = []
with open("z6_e1_1_e2_10.txt") as file:
	for line in file.readlines():
			arr.append([float(i) for i in line.split(' ') if i != '\n'])
fig, ax = plt.subplots()
plt.xlabel('x')
plt.ylabel('y')
plt.title('e1=1,e2=10')
x_max = 100
y_max = 100
plt.xticks(np.arange(0, x_max/10 + 1, 5))
plt.yticks(np.arange(0, y_max/10 + 1, 5))
plt.imshow(arr, cmap='seismic',vmin=-1,vmax=1, extent=[0, x_max/10, 0, y_max/10])
plt.colorbar(ticks=[-0.8,0,0.8])
plt.savefig('z6_e1_1_e2_10_plot.png')