import matplotlib.pyplot as plt

def read_data_file(file_path):
	data = []
	with open(file_path) as file:
		for line in file.readlines():
			data.append([float(i) for i in line.split(' ') if i != '\n'])
	return data

def plot_first():
	vx = read_data_file("out/zad3_vx.txt")
	vy = read_data_file("out/zad3_vy.txt")
	figure1, ((ax1), (ax2)) = plt.subplots(2, 1, figsize=(5, 5))
	ax1.imshow(vx, cmap='jet')
	ax1.set_xlim(0,400)
	ax1.set_ylim(0,90)
	ax1.title.set_text('Vx')
	ax1.set_xlabel('x')
	ax1.set_ylabel('y')

	ax2.imshow(vy, cmap='jet')
	ax2.set_xlim(0,400)
	ax2.set_ylim(0,90)
	ax2.title.set_text('Vy')
	ax2.set_xlabel('x')
	ax2.set_ylabel('y')

	figure1.colorbar(ax1.imshow(vx, cmap='jet', vmin=-5, vmax=45))
	figure1.colorbar(ax2.imshow(vy, cmap='jet', vmin=-20, vmax=20))

	plt.savefig("res/Vx_Vy.png",bbox_inches='tight',transparent=False)

def plot_second():
	c0 = read_data_file("out/zad5_c0.txt")
	c01 = read_data_file("out/zad5_c0.1.txt")
	xsr0 = read_data_file("out/zad5_xsr0.txt")
	xsr01 = read_data_file("out/zad5_xsr0.1.txt")
	t = read_data_file("out/zad3_tt.txt")
	figure1, ((ax1), (ax2)) = plt.subplots(2, 1, figsize=(5, 10))
	ax1.plot(t, c0, label='D=0.0')
	ax1.plot(t, c01, label='D=0.1')
	ax1.set_xlim([0,0.6])
	ax1.set_ylim([0.5,1.05])
	ax1.title.set_text('c(tn)')
	ax1.set_xlabel('tn')
	ax1.set_ylabel('c(tn)')
	ax1.legend(loc='lower left')

	ax2.plot(t, xsr0, label='D=0.0')
	ax2.plot(t, xsr01, label='D=0.1')
	ax2.set_xlim([0,0.6])
	ax2.set_ylim([0,4])
	ax2.title.set_text('xsr(tn)')
	ax2.set_xlabel('tn')
	ax2.set_ylabel('xsr(tn)')
	ax2.legend(loc='lower left')

	plt.savefig("res/C(tn)_Xsr(tn).png",bbox_inches='tight',transparent=False)

import subprocess

def plot_third():
	subprocess.call(f'gnuplot plot3_0.gp', shell=True)
	subprocess.call(f'gnuplot plot3_01.gp', shell=True)

def main():
	plot_first()
	plot_second()
	plot_third()

if __name__ == '__main__':
	main()