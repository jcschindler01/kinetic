

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gammaln
plt.style.use("classic")



fig0 = plt.figure(0, figsize=(4,4))

eps = 0.01
sq = .5-2*eps

outline = plt.axes([eps,eps, 2*(sq+eps), 2*(sq+eps)])

ax0 = plt.axes([eps,.5+eps,sq,sq])
ax1 = plt.axes([.5+eps,.5+eps,sq,sq])
ax2 = plt.axes([eps,eps,sq,sq])
ax3 = plt.axes([.5+eps,eps,sq,sq])
axx = [ax0,ax1,ax2,ax3]



for ax in axx:
	plt.sca(ax)
	plt.xticks([])
	plt.yticks([])
	plt.xlim(0,1)
	plt.ylim(0,1)
	plt.axis("off")

for ax in [outline]:
	plt.sca(ax)
	plt.xticks([])
	plt.yticks([])
	plt.xlim(0,1)
	plt.ylim(0,1)







fig0.savefig("fig5c.pdf")


