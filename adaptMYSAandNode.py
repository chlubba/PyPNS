import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate

diamArray = np.arange(0.1, 20, 0.1)

diamsMYSA = [           0.,     0.25,   0.5,    0.75,   1.,     16.]
MYSALen = np.multiply([ 0.4,    0.55,   0.6,    0.7,    0.8,    1], 3)

MYSASpline = interpolate.splrep(diamsMYSA, MYSALen, k=1)

diamsNode =             [0.,    0.25,   0.5,    0.75,   1.,     16.]
NodeLen = np.multiply(  [0.4,   0.45,   0.55,   0.6,    0.7,    1], 1)

# z1 = np.polyfit(diamsMYSA, MYSALen, 5)
# p1 = np.poly1d(z1)
# plt.plot(diamArray, p1(diamArray))

MYSASpline = interpolate.splrep(diamsMYSA, MYSALen, k=1)
nodeSpline = interpolate.splrep(diamsNode, NodeLen, k=1)

f, (ax1, ax2, ax3) = plt.subplots(3,1, sharex=True)

ax1.plot(diamArray, interpolate.splev(diamArray, MYSASpline) / interpolate.splev(diamArray, nodeSpline), label='ratio MYSA / node length')
ax1.plot(diamArray, interpolate.splev(diamArray, nodeSpline))
ax1.plot(diamArray, interpolate.splev(diamArray, MYSASpline)/3)
ax2.plot(diamArray, interpolate.splev(diamArray, MYSASpline), label='MYSA length')
ax3.plot(diamArray, interpolate.splev(diamArray, nodeSpline), label='node length')

plt.legend()



# f, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2,3)
#
# ax1.scatter(diams, nodeNodeSep, label='Node Separation')
# ax1.set_title('Node Separation')
# z1 = np.polyfit(diams, nodeNodeSep, 1)
# p1 = np.poly1d(z1)
# ax1.plot(diamArray, p1(diamArray))
#
# ax2.scatter(diams, noLamella, label='Number of Lamella')
# ax2.set_title('Number of lamella')
# z2 = np.polyfit(diams, noLamella, 1)
# p2 = np.poly1d(z2)
# ax2.plot(diamArray, p2(diamArray))
#
# ax3.scatter(diams, nodeDiam, label='Node and MYSA Diam')
# ax3.set_title('Node and MYSA Diameter')
# z3 = np.polyfit(diams, nodeDiam, 2)
# p3 = np.poly1d(z3)
# ax3.plot(diamArray, p3(diamArray))
#
# ax4.scatter(diams, FLUTLen, label='FLUT length')
# ax4.set_title('FLUT length')
# z4 = np.polyfit(diams, FLUTLen, 1)
# p4 = np.poly1d(z4)
# ax4.plot(diamArray, p4(diamArray))
#
# ax5.scatter(diams, STINLen, label='STIN length')
# ax5.set_title('STIN length')
# z5 = np.polyfit(diams, STINLen, 1)
# p5 = np.poly1d(z5)
# ax5.plot(diamArray, p5(diamArray))
#
# ax6.scatter(diams, FLUTDiam, label='FLUT and STIN diameter')
# ax6.set_title('FLUT and STIN diameter')
# z6 = np.polyfit(diams, FLUTDiam, 2)
# p6 = np.poly1d(z6)
# ax6.plot(diamArray, p6(diamArray))

plt.legend()
plt.show()