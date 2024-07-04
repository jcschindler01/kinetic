

import numpy as np
import matplotlib.pyplot as plt
import os
plt.style.use("classic")

##
os.chdir("data")

## 
l2e = np.log2(np.e)

## get data
t = np.genfromtxt("Scan.dat",  unpack=True)[0]
StMab = np.genfromtxt("Scan.dat",  unpack=True)[1]
StMa  = np.genfromtxt("ScanA.dat", unpack=True)[1]
SMab  = np.genfromtxt("Smic.dat",  unpack=True)[1]
SMa   = np.genfromtxt("SmicA.dat", unpack=True)[1]
EA = np.genfromtxt("energyA.dat", unpack=True)[1]
EB = np.genfromtxt("energyB.dat", unpack=True)[1]

with open("speciallist.dat", "r") as f:
	Stau = float(f.readline())
	deff = float(f.readline().split("-")[0])
	qmin = float(f.readline())
	m = float(f.readline())

logd = np.log(140**2)

## process data
logt = np.log10(t+1)
StMab *= l2e
StMa  *= l2e
SMab  *= l2e
SMa   *= l2e
Stau  *= l2e
logd  *= l2e

##
plt.figure(figsize=(6,4))
plt.axes([.083,.09,.9,.88])

##
fsize=16
plt.xlabel(r"$t/T$", labelpad=-12, size=fsize)
plt.ylabel(r"$S$" + " " r"$\rm (bits)$", labelpad=-18, size=fsize)
plt.grid(zorder=0)

##
xmin, xmax = (0,3.1)
xticks = np.array([0,10,100,1000,10000])
plt.xticks(np.log10(xticks+1), [r"$0$", r"$10$", r"$10^2$",r"$10^3$",r"$10^4$"], size=fsize)
plt.xlim(xmin, xmax)

##
ymin, ymax = (0,16)
plt.ylim(ymin, ymax)
plt.yticks([ymin,Stau,logd,ymax],[r"$%s$"%ymin,r"$S(\tau)$",r"$\log d$",r"$%d$"%ymax], size=fsize)


##
plt.plot(logt, StMab, 'b-',  marker='.', zorder=119, label=r"$S^{\tau}_{M_A \otimes M_B}$")
plt.plot(logt, StMa , 'c-',  marker='.', zorder=108, label=r"$S^{\tau}_{M_A \otimes \mathbb{1}_B}$")
plt.plot(logt, SMab , 'b--', marker='o', markerfacecolor='w', markeredgecolor='b', markersize=4, zorder=107, label=r"$S_{M_A \otimes M_B}$")
plt.plot(logt, SMa  , 'c--', marker='o', markerfacecolor='w', markeredgecolor='c', markersize=4, zorder=106, label=r"$S_{M_A \otimes \mathbb{1}_B}$")
leg = plt.legend(
	loc='lower right', fontsize=fsize+1, borderpad=.2, labelspacing=.2, 
	borderaxespad=.25, handlelength=2.5, handletextpad=.5
)



##
bbox = dict(fc='.95', ec='w') 
# plt.annotate(text=r"$\log d = 14.3 \;\; {\rm bits}$", xy=(.02,.98), size=fsize-3, xycoords="axes fraction", ha="left", va="top")
# plt.annotate(text=r"$S(\tau) = %.1f \;\; {\rm bits}$"%(Stau), xy=(.02,.92), size=fsize-3, xycoords="axes fraction", ha="left", va="top")


## energy plot
plt.axes([.605,.45,.37,.225])
plt.plot(logt, EA, ':', c='orangered', marker='.', label=r"$\langle H_A \rangle$")
plt.plot(logt, EB, ':', c='dodgerblue', marker='.', label=r"$\langle H_B \rangle$")
plt.xticks(np.log10(xticks+1), ["" for ss in xticks], size=fsize)
plt.yticks([0,-4],[r"$0$",r"$-4$"])
plt.ylabel(r"$E$", labelpad=-15, size=fsize-3)
plt.xlim(xmin,xmax)
plt.ylim(-4,0)
plt.grid()
leg = plt.legend(
	loc='lower right', fontsize=fsize-3, borderpad=0.2, labelspacing=0, 
	borderaxespad=0.2, handlelength=2, handletextpad=0, 
	ncol=2, columnspacing=0.3,
)


os.chdir("..")

##
plt.savefig("002.pdf")

