import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

numberOfPoints = 2
radius = 10

xLimit = 1000
yLimit = 1000
numGrid = 1000

I = 3./numberOfPoints

def distanceOneElectrode(x, y, electrodePosition):

    distance = np.sqrt((x-electrodePosition[0])**2 + (y-electrodePosition[1])**2)

    # voltage = I*resistanceFromInf
    return distance

def voltageAllElectrodes(x, y, electrodePositions, Ra=1., I=1, bipolar=False):

    """

    Args:
        x: x coordinate of interest
        y: y coordinate of interest
        electrodePositions: ...
        Ra: specific resistance in ohm*cm
        I: point current strength in mA
        bipolar: if True, every uneven electrodePosition will be weighted as cathode

    Returns:

        voltage in V

    """
    sumDistance = 0
    for i in range(electrodePositions.shape[0]):
        if bipolar:
            sign = (-1)**i
        else:
            sign = 1
        electrodePosition = electrodePositions[i,:]
        sumDistance += sign/distanceOneElectrode(x, y, electrodePosition)

    constant = Ra/(4*np.pi)*0.01 # Ohm*cm *0.01 = kOhm*um
    resistanceFromInf = constant*sumDistance # kOhm*um/um = kOhm
    voltage = resistanceFromInf*I

    return voltage


# angles = 2*np.pi/numberOfPoints*np.arange(numberOfPoints)
# positions = np.column_stack((np.cos(angles), np.sin(angles))) # np.zeros(numberOfPoints),

numberOfPoints = 20
angles = np.pi/numberOfPoints*np.arange(numberOfPoints)
positions = np.array([]).reshape(0,2)
radius = 500
for angle in angles:
    positions = np.vstack([positions, [np.cos(angle)*radius, np.sin(angle)*radius], [np.cos(angle)*radius, -np.sin(angle)*radius]])
    # positions = np.column_stack((np.cos(angles), np.sin(angles))) # np.zeros(numberOfPoints),

# plt.plot(positions[:,0], positions[:,1])
# plt.show()

# positions = np.column_stack((np.zeros(4), np.array([-500,0,500,0])))


x = np.linspace(-xLimit, xLimit, numGrid)
y = np.linspace(-yLimit, yLimit, numGrid)

xx, yy = np.meshgrid(x, y)#, sparse=True)
# zs = np.array([np.sin(xx**2 + yy**2) / (xx**2 + yy**2) for x,y in zip(np.ravel(xx), np.ravel(yy))])
# zz = zs.reshape(xx.shape)
# z = np.sin(xx**2 + yy**2) / (xx**2 + yy**2)
z = voltageAllElectrodes(xx, yy, positions, Ra=12500, I=I, bipolar=True)
# h = plt.contourf(x,y,z)

limit = 1

zTrunc = np.clip(z, -limit, limit)

ax.plot_surface(xx, yy, zTrunc, cmap=cm.jet)
ax.set_zlim3d(-limit, limit)
plt.show()



