

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gammaln
plt.style.use("classic")


def systemA(x,y,vx,vy,sys=""):
	##
	NN, kk = len(x), np.arange(len(x))
	mask = kk < NN/2
	##
	if sys=="hotcold":
		mask = x<0.5
	if sys=="corner":
		mask = kk <= NN/4
	if sys=="gun":
		mask = x>0.7
	if sys=="chain":
		mask = np.logical_and(kk<NN/4, kk%1 == 0)
	return mask

def sysname(f):
	name = (f.split("/txt/")[-1]).split("_")[0]
	return name

f1 = "../../../runs/txt/hotcold_naive_N500_T5_1720045485.txt"
f2 = "../../../runs/txt/corner_naive_N500_T5_1720110742.txt"
f3 = "../../../runs/txt/gun_naive_N500_T4_1707926341.txt"
f4 = "../../../runs/txt/chain_naive_N500_T5_1720109181.txt"

files = [f1,f2,f3,f4]

t1 = .5

x0  = []
y0  = []
vx0 = []
vy0 = []

x1  = []
y1  = []
vx1 = []
vy1 = []

N = []
maskA = []

for F in files:
	with open(F, 'r') as f:
		##
		f.readline() ## header
		line = f.readline().strip().split()
		data = np.array([float(s.strip().replace(',','')) for s in line[9:]])
		t = data[0]
		k  = data[2::5]
		x  = data[3::5]
		y  = data[4::5]
		vx = data[5::5]
		vy = data[6::5]
		##
		x0  += [x,]
		y0  += [y,]
		vx0 += [vx,]
		vy0 += [vy,]
		N += [len(x),]
		maskA += [systemA(x,y,vx,vy,sys=sysname(F))
]
		##
		while t<t1:
			##
			line = f.readline().strip().split()
			data = np.array([float(s.strip().replace(',','')) for s in line[9:]])
			t = data[0]
			k  = data[2::5]
			x  = data[3::5]
			y  = data[4::5]
			vx = data[5::5]
			vy = data[6::5]
		##
		x1  += [x,]
		y1  += [y,]
		vx1 += [vx,]
		vy1 += [vy,]


fig0 = plt.figure(0, figsize=(4,4))

eps = 0.01
sq = .5-2*eps

ax0 = plt.axes([eps,.5+eps,sq,sq])
ax1 = plt.axes([.5+eps,.5+eps,sq,sq])
ax2 = plt.axes([eps,eps,sq,sq])
ax3 = plt.axes([.5+eps,eps,sq,sq])
axx = [ax0,ax1,ax2,ax3]


fig1 = plt.figure(1, figsize=(4,4))

bx0 = plt.axes([eps,.5+eps,sq,sq])
bx1 = plt.axes([.5+eps,.5+eps,sq,sq])
bx2 = plt.axes([eps,eps,sq,sq])
bx3 = plt.axes([.5+eps,eps,sq,sq])
bxx = [bx0,bx1,bx2,bx3]

plt.figure(0)

for ax in axx+bxx:
	plt.sca(ax)
	plt.xticks([])
	plt.yticks([])
	plt.xlim(0,1)
	plt.ylim(0,1)

sty0 = dict(c='k', marker='.', markersize=2, ls='none',zorder=100)
sty1 = dict(c='0.7', lw=.5, zorder=50)
styA = dict(c='r', marker='.', alpha=.9, markersize=.2, ls='none',zorder=200)
bbox = dict(fc='w', ec='0.9', pad=2, alpha=1)
fsize = 14

for i in range(4):
	plt.sca(axx[i])
	plt.plot(x0[i], y0[i], **sty0)
	mask = maskA[i]
	plt.plot(x0[i][mask], y0[i][mask], **styA)
	s = .001
	v = np.sqrt(np.sum(vx0[i]**2 + vy0[i]**2))/N[i]
	dx, dy = s*vx0[i]/v, s*vy0[i]/v
	for n in range(N[i]):
		plt.plot([x0[i][n],x0[i][n]-dx[n]], [y0[i][n],y0[i][n]-dy[n]], **sty1)
	plt.annotate(text=r"$\rm IC$ $%d$"%(i+1), xy=(.05,.95), ha="left", va="top", bbox=bbox, size=fsize, zorder=500)

for i in range(4):
	plt.sca(bxx[i])
	plt.plot(x1[i], y1[i], **sty0)
	mask = maskA[i]
	plt.plot(x1[i][mask], y1[i][mask], **styA)	
	s = .001
	v = np.sqrt(np.sum(vx1[i]**2 + vy1[i]**2))/N[i]
	dx, dy = s*vx1[i]/v, s*vy1[i]/v
	for n in range(N[i]):
		plt.plot([x1[i][n],x1[i][n]-dx[n]], [y1[i][n],y1[i][n]-dy[n]], **sty1)




fig0.savefig("fig5a.pdf")
fig1.savefig("fig5b.pdf")


