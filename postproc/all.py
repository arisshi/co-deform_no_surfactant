import numpy as np
from numpy import column_stack
from matplotlib import pyplot as plt
from matplotlib.pyplot import *
from math import pi, log
from scipy import integrate, interpolate
import os

# paths
workdir  = os.getcwd()
outdir   = workdir + '/../src/2_dr01_ca01_bo01/'
outdir   = workdir + '/../src/output/'
#outdir   = workdir + '/../output/go/Ca1em03/Bo1ep00/Ma1ep00/dr5em02/dt1em05/'
#outdir   = workdir + '/../output/stop/Ca1em03/Bo1ep00/Ma0ep00/tstop1-2/dr5em02/dt1em05/'

# filenames
hfile = '/evol_h.txt'       
tfile = '/evol_t.txt'       
rfile = '/evol_r.txt'    
#hfile1 = '/aux_h1.txt'
#hfile2 = '/aux_h2.txt'
qfile = '/evol_q.txt'       
pfile = '/evol_p.txt'       
#vfile1 = '/aux_v1.txt'
#vfile2 = '/aux_v2.txt'
#gfile1 = '/aux_gamma1.txt'
#gfile2 = '/aux_gamma2.txt'
efile = '/evol_e.txt'       
ffile = '/evol_f.txt'       
bfile = '/evol_b.txt'       
#pfile = '/aux_pdyn.txt'


# load data
tdata     = np.loadtxt(outdir + tfile, unpack=True,skiprows=0)
rdata     = np.loadtxt(outdir + rfile, unpack=True,skiprows=0)
hdata    = np.loadtxt(outdir + hfile, unpack=True,skiprows=0)
#hdata1    = np.loadtxt(outdir + hfile1, unpack=True,skiprows=0)
#hdata2    = np.loadtxt(outdir + hfile2, unpack=True,skiprows=0)
qdata     = np.loadtxt(outdir + qfile, unpack=True,skiprows=0)
#vdata1    = np.loadtxt(outdir + vfile1, unpack=True,skiprows=0)
#vdata2    = np.loadtxt(outdir + vfile2, unpack=True,skiprows=0)
#gdata1    = np.loadtxt(outdir + gfile1, unpack=True,skiprows=0)
#gdata2    = np.loadtxt(outdir + gfile2, unpack=True,skiprows=0)
edata    = np.loadtxt(outdir + efile, unpack=True,skiprows=0)
fdata    = np.loadtxt(outdir + ffile, unpack=True,skiprows=0)
bdata    = np.loadtxt(outdir + bfile, unpack=True,skiprows=0)
pdata     = np.loadtxt(outdir + pfile, unpack=True,skiprows=0)


# number of space and time points
J1 = len(tdata[:,0])
N1 = len(tdata[0,:])
J  = J1 - 1
N  = N1 - 1

colors = ('b','r','g','m','orange','dimgray','brown','y')

plt.figure(1)
for count in range(0,5):
	if count == 0:
		ydata = hdata;
		ystr = '$h1$'
	if count == 1:
		ydata = bdata;
		ystr = '$h2$'
	if count == 2:
		ydata = edata;
		ystr = '$h_2^{(1)}$'
	if count == 3:
		ydata = qdata;
		ystr  = '$v_r$'
	if count == 4:
		ydata = pdata;
		ystr  = '$P$'
	for i in range(N1):
			itv = 1
			ic = int(np.ceil(i/itv))
			if ( i % itv == 0 and i < 4 ):
				x = rdata[:,i]
				y = ydata[:,i]
				if count == 0:
					colorset = '-'	
					color = colors[ic % len(colors)]
					plt.subplot(221)
					plot(x,y,colorset,label='$\overline{t}$ = ' + str('%.1f' % tdata[0,i]),linewidth=5.5, color = color)
					legend(frameon = False,fontsize=16)
					xlabel('$r$',fontsize=26)
					ylabel(ystr,fontsize=26)
				if count == 1:
					colorset = '--'	
					color = colors[ic % len(colors)]
					temp = plt.subplot(221)
#					plot(x,y,colorset,label='$\overline{t}$ = ' + str('%.1f' % tdata[0,i]),linewidth=5.5, color = color)
					ylabel('$h_i$',fontsize=26)
					xlabel('$r$',fontsize=26)
#	xlim([0,2]);
#    ylim([-2,3]);
				if count == 2:
					colorset = '-'	
					color = colors[ic % len(colors)]
					temp = plt.subplot(222)
					plot(x,y,colorset,label='$\overline{t}$ = ' + str('%.1f' % tdata[0,i]),linewidth=5.5, color = color)
					ylabel(ystr,fontsize=26)
					xlabel('$r$',fontsize=26)
					legend(frameon = False,fontsize=16)
				if count == 3:
					colorset = '-'	
					color = colors[ic % len(colors)]
					temp = plt.subplot(223)
					plot(x,y,colorset,label='$\overline{t}$ = ' + str('%.1f' % tdata[0,i]),linewidth=5.5, color = color)
					ylabel(ystr,fontsize=26)
					xlabel('$r$',fontsize=26)
					legend(frameon = False,fontsize=16)
				if count == 4:
					colorset = '-'	
					color = colors[ic % len(colors)]
					temp = plt.subplot(224)
					plot(x,y,colorset,label='$\overline{t}$ = ' + str('%.1f' % tdata[0,i]),linewidth=5.5, color = color)
					ylabel(ystr,fontsize=26)
					xlabel('$r$',fontsize=26)
					legend(frameon = False,fontsize=16)


	count = count+1

#plt.switch_backend('wxAgg')
#mng = plt.get_current_fig_manager()
#mng.frame.Maximize(True)
#figsize[0] = 12
#figsize[1] = 9
#plt.rcParams["figure.figsize"] = figsize
plt.show()
