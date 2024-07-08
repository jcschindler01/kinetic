
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gammaln, erf
plt.style.use("classic")


l2e = np.log2(np.e)
logfac2 = lambda n: l2e*gammaln(n+1)

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

def STAU(N=500, alpha=1111, d=2):
	return N*d*np.log2(alpha*np.sqrt(np.e)) - logfac2(N)

def numlines(datafile):
	with open(datafile, "rb") as f:
		num_lines = sum(1 for _ in f)
	return num_lines

def D(p,q):
	return np.nansum(p*np.log2(p/q))


def PSTAR(P, Q, dP=.02):
	delta = np.max(np.abs(P-Q))
	lam = np.min([1., 0.5*dP/delta])
	PSTAR = (1-lam)*P + lam*Q
	return PSTAR

def P_spatial(x, y, b=6):
	N = len(x)
	P = np.histogram2d(x, y, bins=np.linspace(0,1,b+1))[0]/N
	return P
	
def Q_spatial(b=6):
	return np.ones((b,b))/b**2

vmax,  dv  = 5, .3
vvmax, dvv = 5, .1
vedges  = np.arange(-vmax, vmax +0.5*dv,  dv)
vvedges = np.arange(0    , vvmax+0.5*dvv, dvv)

def P_velocity(vx, vy, vedges=vedges, eps=1e-9):
	N = len(vx)
	vx, vy = np.clip(vx, -vmax+eps, vmax-eps), np.clip(vy, -vmax+eps, vmax-eps)
	P = np.histogram2d(vx, vy, bins=vedges)[0]/N
	return P

def Q_velocity(vedges=vedges, s=1):
	## get prob in each bin for exp(-v^2/2s^2), s=sqrt(2E/Nd)
	integral = erf(vedges/(np.sqrt(2)*s))
	q = integral[1:]-integral[:-1]
	qq = np.outer(q,q)
	Q = qq/np.sum(qq)
	return Q

def P_speed(v, vedges=vvedges, eps=1e-9):
	N = len(vx)
	v = np.clip(v, 0, vvmax-eps)
	P = np.histogram(v, bins=vedges)[0]/N
	return P

def Q_speed(vedges=vvedges, s=1):
	## get prob in each bin
	integral = np.exp(-(vedges**2/(2*s**2)))
	q = integral[1:]-integral[:-1]
	Q = q/np.sum(q)
	return Q

f1 = "../../../runs/txt/hotcold_naive_N500_T5_1720045485.txt"
f2 = "../../../runs/txt/corner_naive_N500_T5_1720110742.txt"
f3 = "../../../runs/txt/gun_naive_N500_T5_1720111984.txt"
f4 = "../../../runs/txt/chain_naive_N500_T5_1720109181.txt"


free1 = "../../../runs/txt/hotcold_naive_N500_T5_1720045485.txt"
free2 = "../../../runs/txt/corner_naive_N500_T5_1720110742.txt"
free3 = "../../../runs/txt/gun_naive_N500_T5_1720111984.txt"
free4 = "../../../runs/txt/chain_naive_N500_T5_1720109181.txt"




files = [f1,f2,f3,f4]
labels = ['a','b','c','d']



for ic in [0,1,2,3]:

	datafile = files[ic]
	npts = numlines(datafile)-1
	d = 2

	t = np.nan * np.zeros(npts)
	EE  = np.nan * np.zeros(npts)
	S_spatial = np.nan * np.zeros(npts)
	S_velocity = np.nan * np.zeros(npts)
	S_speed = np.nan * np.zeros(npts)
	S_eAB = np.nan * np.zeros(npts)



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
				maskA = systemA(x,y,vx,vy,sys=sysname(datafile))
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
			E = np.sum(e)
			##
			t[n] = tt
			## spatial cg
			P, Q = P_spatial(x,y), Q_spatial()
			Pstar = PSTAR(P,Q)
			S_spatial[n] = Stau - N * D(Pstar,Q)
			## velocity cg
			P, Q = P_velocity(vx,vy), Q_velocity(s=np.sqrt(2*E/(N*d)))
			Pstar = PSTAR(P,Q)
			S_velocity[n] = Stau - N * D(Pstar,Q)
			## speed cg
			v = np.sqrt(vx**2 + vy**2)
			P, Q = P_speed(v), Q_speed(s=np.sqrt(2*E/(N*d)))
			Pstar = PSTAR(P,Q)
			S_speed[n] = Stau - N * D(Pstar,Q)
			## thermodynamic cg
			S_eAB[n] = Stau-(NA*d/2-1)*np.log2(epsA/EAA)-(NB*d/2-1)*np.log2(epsB/(EBB+1e-12))
			##
			n += 1
			initial = False


	##
	plt.figure(figsize=(4,3))
	plt.axes([.12,.12,.85,.83])

	##
	fsize=16
	plt.xlabel(r"$t/T$", labelpad=-12, size=fsize)
	if ic==0:
		plt.ylabel(r"$S/N$" + " " r"$\rm (bits)$", labelpad=-18, size=fsize)
	plt.grid(zorder=0)


	##
	xmin, xmax = (0,3)
	xticks = np.array([0,1,2,3,4,5])
	plt.xticks(xticks, [r"$%s$"%s for s in xticks], size=fsize)
	plt.xlim(xmin, xmax)


	#
	ymin, ymax = (10,15)
	plt.ylim(ymin, ymax)
	plt.yticks([ymin,ymax],[r"$%s$"%ymin,r"$%s$"%ymax], size=fsize)
	if True:
		plt.yticks([ymin,Stau/N,ymax],[r"$%s$"%ymin,r"$S(\tau)$",r"$%s$"%ymax], size=fsize)


	##
	sty = dict(lw=1)

	##
	plt.plot(t, S_spatial/N, 'c-', zorder=125, label=r"$M_{P(\vec{x})}$", **sty)
	plt.plot(t, S_speed/N, 'g-', zorder=123, label=r"$M_{P(v)}$", **sty)
	plt.plot(t, S_velocity/N, 'm-', zorder=124, label=r"$M_{P(\vec{v})}$", **sty)
	plt.plot(t, S_eAB/N, 'b-', zorder=122, label=r"$M_{E_A} \otimes M_{E_B}$", **sty)


	##
	leg = plt.legend(
		loc='lower right', fontsize=fsize, borderpad=.1, labelspacing=.1, 
		borderaxespad=.2, handlelength=2, handletextpad=.3
	)


	##
	plt.savefig("fig6%s.pdf"%(labels[ic]))





















