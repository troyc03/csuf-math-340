# Homework 3: MATH 340
# Name: Troy Chin
# Date: 09/27/25

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d, CubicSpline, PchipInterpolator
import datetime

# Data
x = np.arange(0,7)
y = np.array([15, 18, 21, 20, 17, 14, 11])
u = np.linspace(x.min(), x.max(), 200)

# Interpolants
linear_interp = interp1d(x, y, kind='linear')  
v_linear = linear_interp(u)

coeffs = np.polyfit(x,y,len(x)-1)

v_poly = np.polyval(coeffs, u)

spline_interp = CubicSpline(x,y, bc_type='not-a-knot')
v_spline = spline_interp(u)

pchip_interp = PchipInterpolator(x,y)
v_pchip = pchip_interp(u)

# Plotting
fig, axs = plt.subplots(2, 2, figsize=(10, 8))

axs[0,0].plot(x, y, 'o', u, v_linear, 'g', linewidth=1.5)
axs[0,0].set_title('Linear Interpolant')
axs[0,0].legend(['Data','Linear'])
axs[0,0].grid(True)

axs[0,1].plot(x, y, 'o', u, v_poly, 'r', linewidth=1.5)
axs[0,1].set_title('Polynomial Interpolant')
axs[0,1].legend(['Data','Poly'])
axs[0,1].grid(True)

axs[1,0].plot(x, y, 'o', u, v_spline, 'c', linewidth=1.5)
axs[1,0].set_title('Cubic Spline Interpolant')
axs[1,0].legend(['Data','Spline'])
axs[1,0].grid(True)

axs[1,1].plot(x, y, 'o', u, v_pchip, 'm', linewidth=1.5)
axs[1,1].set_title('PCHIP')
axs[1,1].legend(['Data','PCHIP'])
axs[1,1].grid(True)

plt.suptitle('Comparison of Interpolants')
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.show()

# Question 2

# Define data array
W = np.array([
    [10, 27, 2001,  5, 10,  4,  8],
    [11, 19, 2001,  7,  4,  5, 11],
    [12,  3, 2001,  8, 12,  6,  4],
    [12, 20, 2001, 10, 14,  8,  7],
    [ 1,  9, 2002, 12, 13, 10,  3],
    [ 1, 23, 2002, 14,  8, 12,  0],
    [ 3,  6, 2002, 16, 10, 13, 10]
    ])

# Interpolant data
x = np.array([-1,-0.5,0,0.5,1])
y = np.array([0,-0.5,1,0.5,0])

# Create datetime array
t = [datetime.date(int(row[2]), int(row[0]), int(row[1])) for row in W]

# Calculate weighted means
tom = W[:,3] + W[:, 4]/16
ben = W[:, 5] + W[:, 6]/16

# Plot results
plt.figure(figsize=(8,5))
plt.plot(t, tom, 'o-', label='Tom')
plt.plot(t, ben, 's--', label='Ben')
plt.xlabel('Date')
plt.ylabel('Weights')
plt.title('Tom and Ben\'s Weights Over Time')
plt.legend()
plt.grid(True)
plt.show()

# Linear interpolation
xx = np.linspace(-1,1,200)
f_linear = interp1d(x,y,kind='linear')
yy_pieceln = f_linear(xx)

# Lagrangian Polynomial Interpolation
coeffs = np.polyfit(x,y,len(x)-1)
yy_poly = np.polyval(coeffs, xx)

# Cubic Spline Interpolation
f_spline = CubicSpline(x,y, bc_type='not-a-knot')
yy_spline = f_spline(xx)

# PCHIP Interpolation
f_pchip = PchipInterpolator(x, y)
yy_pchip = f_pchip(xx)

# Plot results
fig, axs = plt.subplots(2, 2, figsize=(10, 8))

axs[0,0].plot(x, y, 'o', xx, yy_pieceln, 'g', linewidth=1.5)
axs[0,0].set_title('Piecewise Linear')
axs[0,0].legend(['Data','piecelin'])
axs[0,0].grid(True)

axs[0,1].plot(x, y, 'o', xx, yy_poly, 'r', linewidth=1.5)
axs[0,1].set_title('Polynomial')
axs[0,1].legend(['Data','polyinterp'])
axs[0,1].grid(True)

axs[1,0].plot(x, y, 'o', xx, yy_spline, 'c', linewidth=1.5)
axs[1,0].set_title('Cubic Spline')
axs[1,0].legend(['Data','splinetx'])
axs[1,0].grid(True)

axs[1,1].plot(x, y, 'o', xx, yy_pchip, 'm', linewidth=1.5)
axs[1,1].set_title('PCHIP')
axs[1,1].legend(['Data','pchiptx'])
axs[1,1].grid(True)

plt.suptitle('Comparison of Interpolants for [-1, 1]')
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.show()

