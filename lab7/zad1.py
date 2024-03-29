import matplotlib.pyplot as plt

def read_data_file(file_path):
	data = []
	with open(file_path) as file:
		for line in file.readlines():
			data.append([float(i) for i in line.split(' ') if i != '\n'])
	return data

def first_plot():
	psi_minus_one = read_data_file("out/psi_-1000.txt")
	u_minus_one = read_data_file("out/u_-1000.txt")
	v_minus_one = read_data_file("out/v_-1000.txt")
	zeta_minus_one = read_data_file("out/zeta_-1000.txt")
	figure1, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(20, 10))

	contour_plot1 = ax1.contour(psi_minus_one, cmap='inferno', levels=30)
	contour_plot1.set_clim(-55, -50)
	ax1.set_xlim(0, 200)
	ax1.set_ylim(0, 90)
	ax1.title.set_text('Q=-1000, psi(x,y)')
	ax1.set_xlabel('x')
	ax1.set_ylabel('y')
	ax1.set_aspect('equal', adjustable='box')
	contour_plot2 = ax2.contour(zeta_minus_one, cmap='inferno', levels=100)
	contour_plot2.set_clim(-200, 350)
	ax2.set_xlim(0,200)
	ax2.set_ylim(0,90)
	ax2.title.set_text('Q=-1000, zeta(x,y)')
	ax2.set_xlabel('x')
	ax2.set_ylabel('y')
	ax2.set_aspect('equal', adjustable='box')
	ax3.imshow(u_minus_one, cmap='jet')
	ax3.set_xlim(0,200)
	ax3.set_ylim(0,90)
	ax3.title.set_text('Q=-1000, u(x,y)')
	ax3.set_xlabel('x')
	ax3.set_ylabel('y')
	ax4.imshow(v_minus_one, cmap='jet')
	ax4.set_xlim(0,200)
	ax4.set_ylim(0,90)
	ax4.title.set_text('Q=-1000, v(x,y)')
	ax4.set_xlabel('x')
	ax4.set_ylabel('y')

	figure1.colorbar(contour_plot1, ax=ax1)
	figure1.colorbar(contour_plot2, ax=ax2)
	figure1.colorbar(ax3.imshow(u_minus_one, cmap='jet', vmin=-2, vmax=16))
	figure1.colorbar(ax4.imshow(v_minus_one, cmap='jet', vmin=-6, vmax=1))
	plt.savefig("res/Q_-1000.png",bbox_inches='tight',transparent=False)

def second_plot():
	psi_minus_four = read_data_file("out/psi_-4000.txt")
	u_minus_four = read_data_file("out/u_-4000.txt")
	v_minus_four = read_data_file("out/v_-4000.txt")
	zeta_minus_four = read_data_file("out/zeta_-4000.txt")
	figure1, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(20, 10))

	contour_plot1 = ax1.contour(psi_minus_four, cmap='inferno', levels=30)
	contour_plot1.set_clim(-218, -202)
	ax1.set_xlim(0,200)
	ax1.set_ylim(0,90)
	ax1.title.set_text('Q=-4000, psi(x,y)')
	ax1.set_xlabel('x')
	ax1.set_ylabel('y')
	ax1.set_aspect('equal', adjustable='box')
	contour_plot2 = ax2.contour(zeta_minus_four, cmap='inferno', levels=100)
	contour_plot2.set_clim(-700, 1100)
	ax2.set_xlim(0,200)
	ax2.set_ylim(0,90)
	ax2.title.set_text('Q=-4000, zeta(x,y)')
	ax2.set_xlabel('x')
	ax2.set_ylabel('y')
	ax2.set_aspect('equal', adjustable='box')
	ax3.imshow(u_minus_four, cmap='jet')
	ax3.set_xlim(0,200)
	ax3.set_ylim(0,90)
	ax3.title.set_text('Q=-4000, u(x,y)')
	ax3.set_xlabel('x')
	ax3.set_ylabel('y')
	ax4.imshow(v_minus_four, cmap='jet')
	ax4.set_xlim(0,200)
	ax4.set_ylim(0,90)
	ax4.title.set_text('Q=-4000, v(x,y)')
	ax4.set_xlabel('x')
	ax4.set_ylabel('y')

	figure1.colorbar(contour_plot1)
	figure1.colorbar(contour_plot2)
	figure1.colorbar(ax3.imshow(u_minus_four, cmap='jet', vmin=-10, vmax=70))
	figure1.colorbar(ax4.imshow(v_minus_four, cmap='jet', vmin=-14, vmax=4))
	plt.savefig("res/Q_-4000.png",bbox_inches='tight',transparent=False)

def third_plot():
	psi_four = read_data_file("out/psi_4000.txt")
	u_four = read_data_file("out/u_4000.txt")
	v_four = read_data_file("out/v_4000.txt")
	zeta_four = read_data_file("out/zeta_4000.txt")
	figure1, ax1 = plt.subplots(1, 1)
	contour_plot = ax1.contour(psi_four, cmap='inferno', levels=30)
	contour_plot.set_clim(202, 218)
	ax1.set_xlim(0,200)
	ax1.set_ylim(0,90)
	ax1.title.set_text('Q=4000, psi(x,y)')
	ax1.set_xlabel('x')
	ax1.set_ylabel('y')
	ax1.set_aspect('equal', adjustable='box')
	figure1.colorbar(contour_plot)
	plt.savefig("res/Q_4000.png",bbox_inches='tight',transparent=False)

def main():
	first_plot()
	second_plot()
	third_plot()

if __name__ == '__main__':
	main()