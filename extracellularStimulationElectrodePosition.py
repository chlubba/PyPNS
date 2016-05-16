import numpy as np
import matplotlib.pyplot as plt

number_contact_points = 8
angles = 2*np.pi/number_contact_points*np.arange(number_contact_points)

print angles

stimCoords = np.column_stack((np.zeros(number_contact_points), np.cos(angles), np.sin(angles)))

print stimCoords

plt.plot(np.cos(angles), np.sin(angles))
plt.show()
