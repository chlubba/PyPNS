import numpy as np
import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')

numberOfPoints = 20
radius = 10

xLimit = 2
yLimit = 2
res = 0.005

I = 8./numberOfPoints

# def voltageOneElectrode(x, y, electrodePosition, Ra=0.3, I=1):
#     constant = Ra/(4*np.pi)*0.01
#     distance = np.sqrt((x-electrodePosition[0])**2 + (y-electrodePosition[1])**2)
#     resistanceFromInf = constant/distance
#     voltage = I*resistanceFromInf
#     return voltage
#
# def voltageAllElectrodes(x, y, electrodePositions):
#     sumVoltage = 0
#     for i in range(electrodePositions.shape[0]):
#         electrodePosition = electrodePositions[i,:]
#         sumVoltage += voltageOneElectrode(x, y, electrodePosition)
#
#     return sumVoltage

def distanceOneElectrode(x, y, electrodePosition):

    distance = np.sqrt((x-electrodePosition[0])**2 + (y-electrodePosition[1])**2)

    # voltage = I*resistanceFromInf
    return distance

def voltageAllElectrodes(x, y, electrodePositions, Ra=0.3, I=1):
    sumDistance = 0
    for i in range(electrodePositions.shape[0]):
        electrodePosition = electrodePositions[i,:]
        sumDistance += 1/distanceOneElectrode(x, y, electrodePosition)

    constant = Ra/(4*np.pi)*0.01
    resistanceFromInf = constant*sumDistance
    voltage = resistanceFromInf*I

    return voltage


anglesElectrodes = 2*np.pi/numberOfPoints*np.arange(numberOfPoints)
positions = np.column_stack((radius*np.cos(anglesElectrodes), radius*np.sin(anglesElectrodes))) # np.zeros(numberOfPoints),

numberOfPlottingAngles = 400
plottingAngles = 2*np.pi/numberOfPlottingAngles*np.arange(numberOfPlottingAngles)

radiusMax = 8.
radii = np.linspace(0.1, radiusMax, 500)

minVoltage = []
maxVoltage = []

# f, axarr = plt.subplots(2,2)

# for numberOfPoints in range(2,18,4):


for plottingRadius in radii:
    x, y = plottingRadius*np.cos(plottingAngles), plottingRadius*np.sin(plottingAngles)
    voltages = voltageAllElectrodes(x, y, positions, I=I)
    minVoltage.append(min(voltages))
    maxVoltage.append(max(voltages))

    # axarr[0,1].plot(x,y)


# plt.plot(radii, minVoltage, label='min')
# plt.plot(radii, maxVoltage, label='max')
# plt.legend()

# calculate maximum difference in % between min and max and r=0 and r=...

numberOfPointsArray = range(2,40,1) # [2, 40] #

overallMinVoltage = []
overallMaxVoltage = []
radiusMaxMinVoltage = []

for numberOfPoints in numberOfPointsArray:

    anglesElectrodes = 2*np.pi/numberOfPoints*np.arange(numberOfPoints)
    positions = np.column_stack((radius*np.cos(anglesElectrodes), radius*np.sin(anglesElectrodes))) # np.zeros(numberOfPoints),

    minVoltageRadiuswise = []
    maxVoltageRadiuswise = []

    for plottingRadius in radii:
        x, y = plottingRadius*np.cos(plottingAngles), plottingRadius*np.sin(plottingAngles)
        voltages = voltageAllElectrodes(x, y, positions, I=I)
        minVoltageRadiuswise.append(min(voltages))
        maxVoltageRadiuswise.append(max(voltages))

    # plt.plot(radii, minVoltageRadiuswise)
    # plt.plot(radii, maxVoltageRadiuswise)
    # plt.title(str(numberOfPoints)+' points')
    # plt.show()

    overallMinVoltage.append(min(minVoltageRadiuswise))
    overallMaxVoltage.append(max(maxVoltageRadiuswise))
    radiusMaxMinVoltage.append(minVoltageRadiuswise[-1])

# plt.plot(numberOfPointsArray, overallMaxVoltage)
# plt.plot(numberOfPointsArray, overallMinVoltage)

# plt.plot(numberOfPointsArray, np.divide(overallMaxVoltage, overallMinVoltage))
plt.plot(numberOfPointsArray, np.add(np.divide(overallMaxVoltage, radiusMaxMinVoltage),-1))


plt.show()


# x = np.arange(-xLimit, xLimit, res)
# y = np.arange(-yLimit, yLimit, res)
#
# xx, yy = np.meshgrid(x, y)#, sparse=True)
# # zs = np.array([np.sin(xx**2 + yy**2) / (xx**2 + yy**2) for x,y in zip(np.ravel(xx), np.ravel(yy))])
# # zz = zs.reshape(xx.shape)
# # z = np.sin(xx**2 + yy**2) / (xx**2 + yy**2)
# z = voltageAllElectrodes(xx, yy, positions, I=I)
# # h = plt.contourf(x,y,z)
# ax.plot_surface(xx, yy, z)
# ax.set_zlim3d(0, 0.005)
# plt.show()



