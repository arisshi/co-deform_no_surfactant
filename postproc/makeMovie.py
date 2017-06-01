import numpy as np
from numpy import column_stack
#import matplotlib
#matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib import animation as manimation
from math import pi, log
from scipy import integrate, interpolate
import os

# paths
workdir  = os.getcwd()
outdir   = workdir + '/../src/output/'
outname = outdir + 'ca1.mp4'
Tstop    = 12
#outdir   = workdir + '/../output/go/Ca1em01/Bo1ep00/Ma1ep01/dr5em02/dt1em05/'
#outdir   = workdir + '/../output/go/Ca1em02/Bo1ep00/Ma0ep00/dr5em02/dt1em05/'
#outdir   = workdir + '/../output/stop/Ca1em03/Bo1ep00/Ma0ep00/tstop1-2/dr5em02/dt1em04/'
#outdir   = workdir + '/../output/stop/Ca1em01/Bo1ep00/Ma1ep02/tstop1-2/dr5em02/dt1em04/'

# filenames
tfile    = '/evol_t.txt'
rfile    = '/evol_r.txt'
hfile    = '/evol_h.txt'
efile    = '/evol_e.txt'
qfile    = '/evol_q.txt'
pfile    = '/evol_p.txt'
#gfile    = '/evol_g.txt'
ffile    = '/evol_f.txt'
bfile    = '/evol_b.txt'
#pfile    = '/evol_p.txt'
#vfile    = '/evol_vs.txt'

# load data
tdatain  = np.loadtxt(outdir + tfile, unpack=True,skiprows=0)
rdatain  = np.loadtxt(outdir + rfile, unpack=True,skiprows=0)
hdatain  = np.loadtxt(outdir + hfile, unpack=True,skiprows=0)
edatain  = np.loadtxt(outdir + efile, unpack=True,skiprows=0)
qdatain  = np.loadtxt(outdir + qfile, unpack=True,skiprows=0)
pdatain  = np.loadtxt(outdir + pfile, unpack=True,skiprows=0)
fdatain  = np.loadtxt(outdir + ffile, unpack=True,skiprows=0)
bdatain  = np.loadtxt(outdir + bfile, unpack=True,skiprows=0)
#gdata    = np.loadtxt(outdir + gfile, unpack=True,skiprows=0)
#pdata    = np.loadtxt(outdir + pfile, unpack=True,skiprows=0)
#vdata    = np.loadtxt(outdir + vfile, unpack=True,skiprows=0)

tdata    = tdatain[:,0:Tstop]
rdata    = rdatain[:,0:Tstop]
hdata    = hdatain[:,0:Tstop]
edata    = edatain[:,0:Tstop]
fdata    = fdatain[:,0:Tstop]
bdata    = bdatain[:,0:Tstop]
qdata    = qdatain[:,0:Tstop]
pdata    = pdatain[:,0:Tstop]

J1 = len(tdata[:,0])
N1 = len(tdatain[0,:])

# set up the figure, the axis, and the plot element to be animated
fig = plt.figure(figsize=(15,4))
ax1 = fig.add_subplot(131)
ax2 = fig.add_subplot(132)
ax3 = fig.add_subplot(133)
line1a, = ax1.plot([], [], 'b-')
line1b, = ax1.plot([], [], 'k-')
line2,  = ax2.plot([], [], 'k-')
line3,  = ax3.plot([], [], 'k-')
time1    = ax1.text(0.8, 0.8, '', transform=ax1.transAxes)
time2    = ax2.text(0.8, 0.8, '', transform=ax2.transAxes)
time3    = ax3.text(0.8, 0.8, '', transform=ax3.transAxes)

for ax in (ax1,ax2,ax3):
#	tickes = np.arange(0,2.1,0.5)
#	ax.set_xticks(tickes)
	ax.set_xlim(0, 10)
	ax.set_xlabel('$\overline{r}$')

ax1.set_ylabel('$\overline{h}$')
ax1.set_ylim(-10, 2)
#tickes = np.arange(-2,2.1,1)
#ax1.set_yticks(tickes)
ax2.set_ylabel('$\overline{v_r}$')
ax2.set_ylim( 0, 1.2)
#tickes = np.arange(0,0.21,0.05)
#ax2.set_yticks(tickes)
ax3.set_ylabel('$\overline{P}$')
ax3.set_ylim( -1, 5)
#tickes = np.arange(-1,1.1,0.5)
#ax3.set_yticks(tickes)


# initialization function: plot the background of each frame
def init() :
	line1a.set_data([], [])
	line1b.set_data([], [])
	line2 .set_data([], [])
	line3 .set_data([], [])
	time1 .set_text('')
	time2 .set_text('')
	time3 .set_text('')
	return line1a, line1b, line2, line3, time1, time2, time3

# animation function, to be called sequentially
def animate(i) :
	t   =  tdata[0,i]
	x   =  rdata[:,i]
	y1a =  hdata[:,i]
	y1b =  bdata[:,i]
	y2  =  qdata[:,i]
	y3  =  pdata[:,i]

	line1a.set_data(x, y1a)
	line1b.set_data(x, y1b)
	line2 .set_data(x, y2 )
	line3 .set_data(x, y3 )
	time1 .set_text("t = " + str(t))
	time2 .set_text("t = " + str(t))
	time3 .set_text("t = " + str(t))
	return line1a, line1b, line2, line3, time1, time2, time3


# call the animator. blit=True means only re-draw the parts that have changed
anim = manimation.FuncAnimation(fig, animate, init_func=init, frames=Tstop, interval=500, repeat=False, blit=True)

# format and save movie
#writer = manimation.FFMpegWriter(fps=5)
#anim.save(outname, writer=writer)

plt.show()
