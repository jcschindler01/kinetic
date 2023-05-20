
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as an
#plt.style.use("classic")


###### read data ########

## alert
print("..reading")

data = np.genfromtxt("data.txt", delimiter=",", skip_header=1, unpack=True)

#i = data[0].astype(int)
t = data[6].astype(float)
i = np.arange(len(t)).astype(int)

x = data[9:-1:5].astype(float)
y = data[10:-1:5].astype(float)
vx = data[11:-1:5].astype(float)
vy = data[12:-1:5].astype(float)

n = len(i)
N = len(x)
dt = t[1]-t[0]

v = np.sqrt(vx**2 + vy**2)
E = 0.5*np.sum(v**2,axis=0)

vmax = 1.0*np.ceil(np.max(v))

########################


########## animate ##########
if True:

	## alert
	print("..animating")

	## frames
	frames = np.arange(len(i))[::1]

	## style
	sty = dict(markersize=1)

	## fig
	fig, ax, = plt.subplots(2,2, layout="constrained", figsize=(5,5))
	ax = ax.flatten()

	## box
	plt.sca(ax[0])
	ax[0].set_aspect(1)
	ax[0].set_xlim(0,1)
	ax[0].set_ylim(0,1)
	ax[0].set_xticks([0,1])
	ax[0].set_yticks([0,1])
	ax[0].set_xlabel("$x$")
	ax[0].set_ylabel("$y$")

	## velocities
	ax[1].set_aspect(1)
	ax[1].set_xlim(-vmax,vmax)
	ax[1].set_ylim(-vmax,vmax)
	ax[1].set_xticks(np.arange(-vmax,vmax+1))
	ax[1].set_yticks(np.arange(-vmax,vmax+1))
	ax[1].set_xlabel("$v_x$")
	ax[1].set_ylabel("$v_y$")

	## histogram
	#ax[2].set_aspect(1)
	ax[2].set_xlim(0,2)
	ax[2].set_ylim(0,2) 
	ax[2].set_xticks(np.arange(vmax+1))
	ax[2].set_yticks([0,1,2])
	ax[2].set_xlabel("$v$")
	ax[2].set_ylabel("$p(v)$")
	nbins = 75

	## others
	for axx in [ax[3]]:
		axx.set_aspect(1)
		axx.set_xlim(0,1)
		axx.set_ylim(0,1)
		axx.set_xticks([])
		axx.set_yticks([])

	## annotation
	note =  ax[3].annotate(text="t= \nE= ", xy=[.5,.5], ha="center", va="bottom")

	## init artists
	xdata,  ydata  = [], []
	vxdata, vydata = [], []
	vdata = []
	particles,  = ax[0].plot([], [], 'ko', **sty)
	velocities, = ax[1].plot([], [], 'ko', **sty)
	bins = np.linspace(0,vmax,nbins)
	bincenter = bins[0:-1] + 0.5*bins[1]
	hist, = ax[2].plot([],[], 'g*')

	## pass artists
	def init():
	    return particles, velocities, hist, note,

	## update artists
	def update(frame):
	    xdata,  ydata  =  x[:, frame],  y[:, frame]
	    particles.set_data(xdata, ydata)
	    vxdata, vydata = vx[:, frame], vy[:, frame]
	    velocities.set_data(vxdata, vydata)
	    vdata, bins0 = np.histogram(v[:, frame], bins=bins, density=True)
	    hist.set_data(bincenter, vdata)
	    note.set_text("T = %.3f\n\nE/N = %.3f\n\nN = %d"%(t[frame],E[frame]/N,N))
	    return particles, velocities, hist, note,

	## animation
	ani = an.FuncAnimation(
			fig, update, init_func=init, blit=True,
			frames=frames, 
			interval=50, 
		)

	## save
	ani.save("out.gif", dpi=200)

	## show
	#plt.show()


#############################





