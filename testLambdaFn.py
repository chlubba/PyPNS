import numpy as np

peakFactor = lambda angle, xP: max(0, (1 - np.mod(angle, np.pi) / np.pi * 5)) * np.min([1, (xP / 0.000190) ** 5])

a = 1.9E-9  # 2.5E-9;
b = 0.00005

# zValues = np.arange(0, 1.05, 0.01) / 100
vValues = lambda zValues, angle, xP: a * (1.0 / (np.abs(zValues) + b)) * peakFactor(angle, xP) + np.maximum(0, (
np.subtract(1, np.abs(zValues / 0.01)) * 8.83e-5))

print vValues(0.001,0,0)