
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gammaln
plt.style.use("classic")

## 
datafile = "../../../runs/txt/hotcold_naive_N500_T1_1711556349.txt"

##
l2e = np.log2(np.e)
logfac2 = lambda n: l2e*gammaln(n+1)


##
def systemA(x,y,vx,vy):
	mask = x < 0.5
	return mask

##
def STAU(N=500, alpha=1111, d=2):
	return N*d*np.log2(alpha*np.sqrt(np.e)) - logfac2(N)

##
def numlines(file):
	with open(datafile, "rb") as f:
		num_lines = sum(1 for _ in f)
	return num_lines

##
npts = numlines(datafile)-1
t = np.nan * np.zeros(npts)
StMab = np.nan * np.zeros(npts)
StMa  = np.nan * np.zeros(npts)
SMab  = np.nan * np.zeros(npts)
SMa   = np.nan * np.zeros(npts)
EE  = np.nan * np.zeros(npts)
EA = np.nan * np.zeros(npts)
EB = np.nan * np.zeros(npts)
d = 2


##
with open(datafile, 'r') as f:
	## header (t=heads[9], x=heads[12,17], (k x y vx vy))
	heads = f.readline().strip().split(',')
	initial = True
	n = 0
	##
	for line in f:
		## get data
		line = line.strip().split()
		data = np.array([float(s.strip().replace(',','')) for s in line[9:]])
		tt = data[0]
		k  = data[2::5]
		x  = data[3::5]
		y  = data[4::5]
		vx = data[5::5]
		vy = data[6::5]
		e = 0.5*(vx**2 + vy**2)
		##
		if initial:
			##
			maskA = systemA(x,y,vx,vy)
			maskB = ~maskA
			##
			E0 = np.sum(e)
			N = len(k)
			NA, NB = np.sum(maskA), np.sum(maskB)
			kA, kB = k[maskA], k[maskB]
			epsA, epsB = NA*E0/N, NB*E0/N
			Stau = STAU(N)
			beta = N*d/(2*E0)
			##
			x0, y0, vx0, vy0 = 1*x, 1*y, 1*vx, 1*vy
		##
		eA, eB = e[maskA], e[maskB]
		EAA, EBB = np.sum(eA), np.sum(eB)
		##
		t[n] = tt
		StMab[n] = Stau-(NA*d/2-1)*np.log2(epsA/EAA)-(NB*d/2-1)*np.log2(epsB/EBB)
		StMa[n] = Stau-(NA*d/2-1)*np.log2(epsA/EAA)-beta*(EAA-epsA)*l2e
		SMab[n] = np.nan
		SMa[n]  = Stau-(NA*d/2-1)*np.log2(epsA/EAA)
		EE[n] = np.sum(e)
		EA[n] = EAA
		EB[n] = EBB
		##
		n += 1
		initial = False


##
plt.figure(figsize=(6,4))
plt.axes([.083,.09,.9,.88])

##
fsize=16
plt.xlabel(r"$t/T$", labelpad=-12, size=fsize)
plt.ylabel(r"$S/N$" + " " r"$\rm (bits)$", labelpad=-18, size=fsize)
plt.grid(zorder=0)


##
xmin, xmax = (0,3.1)
xticks = np.array([0,1,2,3])
plt.xticks(xticks, [r"$%s$"%s for s in xticks], size=fsize)
plt.xlim(xmin, xmax)


##
ymin, ymax = (11,15)
plt.ylim(ymin, ymax)
plt.yticks([ymin,Stau/N,ymax],[r"$%s$"%ymin,r"$S(\tau)$",r"$%d$"%ymax], size=fsize)


##
plt.plot(t, StMab/N, 'b-',  lw=2, zorder=119, label=r"$S^{\tau}_{M_A \otimes M_B}$")
plt.plot(t, StMa/N , 'c-',  lw=2, zorder=108, label=r"$S^{\tau}_{M_A \otimes \mathbb{1}_B}$")
plt.plot(t, SMab/N , 'b-',  lw=2, zorder=107, label=r"$S_{M_A \otimes M_B}$")
plt.plot(t, SMa/N  , 'c--', lw=2, zorder=106, label=r"$S_{M_A} + c$")
leg = plt.legend(
	loc='lower right', fontsize=fsize+1, borderpad=.2, labelspacing=.2, 
	borderaxespad=.25, handlelength=2.5, handletextpad=.5
)



##
bbox = dict(fc='.95', ec='w') 

## energy plot
plt.axes([.605,.45,.37,.225])
plt.plot(t, EA/N, '-', c='orangered', label=r"$E_A$")
plt.plot(t, EB/N, '-', c='dodgerblue', label=r"$E_B$")
plt.xticks(xticks, ["" for ss in xticks], size=fsize)
plt.yticks([0,.5],[r"$0$", r"$.5$"])
plt.ylabel(r"$E/N$", labelpad=-15, size=fsize-3)
plt.xlim(xmin,xmax)
plt.ylim(0,.5)
plt.grid()
leg = plt.legend(
	loc='lower right', fontsize=fsize-3, borderpad=0.1, labelspacing=0, 
	borderaxespad=0.2, handlelength=1.5, handletextpad=0, 
	ncol=2, columnspacing=0.15,
)


## ic plot
faspect = 6./4
width = .25
plt.axes([.25,.2, width,faspect*width])
plt.xticks([])
plt.yticks([])
plt.xlim(0,1)
plt.ylim(0,1)
plt.plot(x0[maskA], y0[maskA], 'r.', markersize=2)
plt.plot(x0[maskB], y0[maskB], 'b.', markersize=2)
bbox = dict(ec='w', fc='w', pad=1, alpha=.9)
plt.annotate(text=r"$A$", xy=(.25, -.03), xycoords="axes fraction", ha='center', va='top', size=fsize-1)
plt.annotate(text=r"$B$", xy=(.75, -.03), xycoords="axes fraction", ha='center', va='top', size=fsize-1)
plt.annotate(text=r"$\rm Initial$", xy=(.5, 1.03), xycoords="axes fraction", ha='center', va='bottom', size=fsize-2)
plt.annotate(text=r"$\rm hot$", xy=(.25, .5), xycoords="axes fraction", ha='center', va='center', rotation=0, bbox=bbox, size=fsize-2)
plt.annotate(text=r"$\rm cold$", xy=(.75, .5), xycoords="axes fraction", ha='center', va='center', rotation=-0, bbox=bbox, size=fsize-2)


##
plt.savefig("fig4.pdf")

