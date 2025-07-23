"""
Symplectic (semi-implicit) Euler for the forced oscillator:
y''(t) + y(t) = sin(t),   y(0) = 0,  y'(0) = 1
Integrated as the first-order system:
    x' = v
    v' = -x + sin(t)
Plots position and velocity vs time.
"""

import numpy as np
import matplotlib.pyplot as plt

# Parameters
T = 10.0
h = 0.01
N = int(T / h) + 1
t = np.linspace(0, T, N)

# Initial conditions
x = np.zeros(N)
v = np.zeros(N)
x[0] = 0
v[0] = 1

# Symplectic Euler integration
for n in range(1, N):
    x[n] = x[n-1] + h * v[n-1]
    v[n] = v[n-1] + h * (-x[n-1] + np.sin(t[n-1]))

# Plot position and velocity
plt.figure(figsize=(10,4))
plt.subplot(1,2,1)
plt.plot(t, x, label='x (position)')
plt.xlabel('Time t')
plt.ylabel('Position x')
plt.grid(True)
plt.title('Position vs Time')

plt.subplot(1,2,2)
plt.plot(t, v, color='orange', label='v (velocity)')
plt.xlabel('Time t')
plt.ylabel('Velocity v')
plt.grid(True)
plt.title('Velocity vs Time')

plt.tight_layout()
plt.show()