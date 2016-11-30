import matplotlib.pyplot as plt
import numpy as np

diams = [5.7, 7.3, 8.7, 10., 11.5, 12.8, 14., 15., 16.]
nodeNodeSep = [500, 760, 1000, 1150, 1250, 1350, 1400, 1450, 1500]
noLamella = [80, 100, 110, 120, 130, 135, 140, 145, 150]
nodeDiam = [1.9, 2.4, 2.8, 3.3, 3.7, 4.2, 4.7, 5., 5.5]
MYSADiam = nodeDiam
FLUTLen = [35, 38, 40, 46, 50, 54, 56, 58, 60]
FLUTDiam = [3.4, 4.6, 5.8, 6.9, 8.1, 9.2, 10.4, 11.5, 12.7]
STINLen = [70.5, 111.2, 152.2, 175.2, 190.5, 205.8, 213.5, 221.2, 228.8]
STIDiam = FLUTDiam

diamArray = np.arange(0.1, 20, 0.1)

f, ((ax1, ax2, ax3)) = plt.subplots(1,3, sharex=True)#, figsize=(20,10))

ax2.scatter(diams, nodeNodeSep, label='Node Separation')
ax2.set_title('Node separation')
z1 = np.polyfit(diams, nodeNodeSep, 1)
p1 = np.poly1d(z1)
ax2.plot(diamArray, p1(diamArray))


ax3.scatter(diams, noLamella, label='Number of Lamella')
ax3.set_title('Number of lamella')
z2 = np.polyfit(diams, noLamella, 1)
p2 = np.poly1d(z2)
ax3.plot(diamArray, p2(diamArray))

ax1.scatter(diams, nodeDiam)
ax1.set_title('Diameter')
z3 = np.polyfit(diams, nodeDiam, 2)
p3 = np.poly1d(z3)
ax1.plot(diamArray, p3(diamArray), label='Node and MYSA')

ax1.scatter(diams, FLUTDiam)
# ax1.set_title('FLUT and STIN diameter')
z6 = np.polyfit(diams, FLUTDiam, 2)
p6 = np.poly1d(z6)
ax1.plot(diamArray, p6(diamArray), label='FLUT and STIN')
ax1.legend(loc='best')
ax1.set_ylim(ymin=0)

ax3.set_xlim([0,20])


# plt.legend()
plt.tight_layout()

import matplotlib2tikz as mtz
mtz.save('extrapolateMcIntyre3.tex', figureheight = '\\figureheight',
           figurewidth = '\\figurewidth')

plt.show()

