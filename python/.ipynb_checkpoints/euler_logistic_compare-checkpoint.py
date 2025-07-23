"""
Explicit Euler for the Logistic Growth IVP
y'(t) = (1 - y/100) y,   y(0) = 1,   t ∈ [0, 10]
Compares numerical solution with exact solution for various step sizes.
"""

import numpy as np
import matplotlib.pyplot as plt

def logistic_rhs(y):
    return (1 - y/100) * y

def exact_solution(t):
    return 100 / (1 + 99 * np.exp(-t))

def euler_logistic(y0, t):
    y = np.zeros_like(t)
    y[0] = y0
    h = t[1] - t[0]
    for n in range(1, len(t)):
        y[n] = y[n-1] + h * logistic_rhs(y[n-1])
    return y

tvals = [0.2, 0.1, 0.05, 0.01]
y0 = 1
T = 10

plt.figure(figsize=(9,6))

for h in tvals:
    t = np.arange(0, T+h, h)
    y_num = euler_logistic(y0, t)
    y_ex = exact_solution(t)
    plt.plot(t, y_num, '--', label=f'Explicit Euler (h={h})')
    if h == min(tvals):  # Plot exact only once
        plt.plot(t, y_ex, 'k-', label='Exact Solution', linewidth=2)

plt.xlabel('t')
plt.ylabel('y(t)')
plt.title('Logistic Growth IVP: Explicit Euler vs. Exact')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()