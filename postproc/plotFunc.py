import numpy as np
from numpy import column_stack
from matplotlib import pyplot as plt 
from matplotlib.pyplot import *
from math import pi, log
from scipy import integrate, interpolate
import os

# paths
workdir  = os.getcwd()
outdir   = workdir + '/../src/output'
#outdir   = workdir + '/../output/go/Ca1em03/Bo1ep00/Ma1ep00/dr5em02/dt1em05/'
#outdir   = workdir + '/../output/stop/Ca1em03/Bo1ep00/Ma0ep00/tstop1-2/dr5em02/dt1em05/'

# filenames
hfile = '/evol_h.txt'       
#hbfile = '/evol_hb.txt'       
qfile = '/evol_q.txt'       
#'/evol_f.txt'       
tfile = '/evol_t.txt'       
rfile = '/evol_r.txt'    
#hfile1 = '/aux_h1.txt'
#hfile2 = '/aux_h2.txt'
#vfile1 = '/aux_v1.txt'
#vfile2 = '/aux_v2.txt'
efile = '/evol_e.txt'       
pfile = '/evol_p.txt'       
#pfile = '/aux_pdyn.txt'
#gfile1 = '/aux_gamma1.txt'
#gfile2 = '/aux_gamma2.txt'


# load data
tdata     = np.loadtxt(outdir + tfile, unpack=True,skiprows=0)
rdata     = np.loadtxt(outdir + rfile, unpack=True,skiprows=0)
hdata    = np.loadtxt(outdir + hfile, unpack=True,skiprows=0)
#hbdata    = np.loadtxt(outdir + hbfile, unpack=True,skiprows=0)
qdata     = np.loadtxt(outdir + qfile, unpack=True,skiprows=0)
#hdata1    = np.loadtxt(outdir + hfile1, unpack=True,skiprows=0)
#hdata2    = np.loadtxt(outdir + hfile2, unpack=True,skiprows=0)
#vdata1    = np.loadtxt(outdir + vfile1, unpack=True,skiprows=0)
#vdata2    = np.loadtxt(outdir + vfile2, unpack=True,skiprows=0)
edata    = np.loadtxt(outdir + efile, unpack=True,skiprows=0)
pdata     = np.loadtxt(outdir + pfile, unpack=True,skiprows=0)
#gdata1    = np.loadtxt(outdir + gfile1, unpack=True,skiprows=0)
#gdata2    = np.loadtxt(outdir + gfile2, unpack=True,skiprows=0)

pl = 'r'
ystr = 'h'

if ystr == 'gap':
	ydata = hdata - edata
	ystr = 'Film Thickness = h_1 - h_2'
if ystr =='e':
	ydata = edata;
	ystr  = '$\overline{\Gamma}_1+\overline{\Gamma}_2$'
if ystr =='hb':
	ydata = hbdata;
	ystr  = '$\overline{h}_1-\overline{h}_2$'
if ystr =='h':
	ydata = hdata/0.0001;
	ystr  = '$\overline{h}_1+\overline{h}_2$'
if ystr =='h1':
	ydata = (hdata + edata )/2.0;
	ystr  = '$\overline{h}_1$'
if ystr =='h2':
	ydata = -(hdata - edata )/2.0;
	ystr  = '$\overline{h}_2$'
if ystr =='p':
	ydata = pdata
	ystr  = '$\overline{p}$'
if ystr =='q':
	ydata = qdata
	ystr  = '$\overline{q}$'
if ystr =='v1':
	ydata = vdata1
	ystr  = '$\overline{v}_{1s}$'
if ystr =='v2':
	ydata = vdata2
	ystr  = '$\overline{v}_{2s}$'
if ystr =='g1':
	ydata = gdata1
	ystr  = '$\overline{\Gamma}_{1}$'
if ystr =='g2':
	ydata = gdata2
	ystr  = '$\overline{\Gamma}_{2}$'

# number of space and time points
J1 = len(tdata[:,0])
N1 = len(tdata[0,:])
J  = J1 - 1
N  = N1 - 1

if pl == 'r' :
	for i in range(N1):
		if ( i % 3 == 0 ):
			x = rdata[:,i]
			y = ydata[:,i]

			plot(x,y,'-',label='$\overline{t}$ = ' + str(tdata[0,i]))


	legend(frameon=False)
	xlabel('$\overline{\sigma}$')
	ylabel(ystr)


if pl == 't':
	for i in range(J1) :
		if i == 0 : # centerline
			x = tdata[i,:]
			y = ydata[i,:]
			plot(x,y,'-')
	xlabel('$\overline{t}$')
	ylabel(ystr)
show()
	


