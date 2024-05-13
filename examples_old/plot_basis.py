import numpy as np
import matplotlib.pyplot as plt
plt.style.use('fivethirtyeight')
basis_p = np.genfromtxt('examples/basis_p.csv', delimiter=',')
basis_v = np.genfromtxt('examples/basis_v.csv', delimiter=',')
basis_a = np.genfromtxt('examples/basis_a.csv', delimiter=',')

plt.figure(figsize=(10, 5))
t = np.linspace(0, 1, len(basis_p[0]))
for y in basis_p:
    plt.plot(t, y, solid_capstyle='round')

plt.figure(figsize=(10, 5))
for y in basis_v:
    plt.plot(t, y, solid_capstyle='round')

plt.figure(figsize=(10, 5))
for y in basis_a:
    plt.plot(t, y, solid_capstyle='round')
plt.show()