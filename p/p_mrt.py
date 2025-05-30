import numpy as np
import matplotlib.pyplot as plt


data = np.loadtxt("output/inte_r/0/150.txt", skiprows=1)
m = data[:, 0]
r = data[:, 1]
y = data[:, 2]

GMsun_c2km = 1.476625038507424

c = m / r * GMsun_c2km

if False:
    Lamb = (
        16
        / 15
        * (1 - 2 * c) ** 2
        * (2 * c * (y - 1) - y + 2)
        / (
            2
            * c
            * (
                4 * (y + 1) * c**4
                + (6 * y - 4) * c**3
                + (26 - 22 * y) * c**2
                + 3 * (5 * y - 8) * c
                - 3 * y
                + 6
            )
            + 3 * (1 - 2 * c) ** 2 * (2 * c * (y - 1) - y + 2) * np.log(1 - 2 * c)
        )
    )
else:
    Lamb = (
        16
        / 15
        * (1 - 2 * c) ** 2
        * (2 - 2 * c + (2 * c - 1) * y)
        / (
            (3 - 12 * c + 13 * c**2 - 2 * c**3 + 2 * c**4) * 4 * c
            + (-3 + 15 * c - 22 * c**2 + 6 * c**3 + 4 * c**4) * 2 * c * y
            + 3 * (1 - 2 * c) ** 2 * (2 - 2 * c + (2 * c - 1) * y) * np.log(1 - 2 * c)
        )
    )

plt.figure(figsize=(12, 4))
plt.subplot(1, 3, 1)
plt.plot(r, m, ".", ms=1)
plt.xlabel("Radius (km)")
plt.ylabel("Mass (M$\\odot$)")

plt.subplot(1, 3, 2)
plt.plot(r, Lamb, ".", ms=1)
plt.xlabel("Radius (km)")
plt.ylabel(r"$\Lambda$")
plt.ylim(0, 1000)

plt.subplot(1, 3, 3)
plt.plot(m, Lamb, ".", ms=1)
plt.xlabel("Mass (M$\\odot$)")
plt.ylabel(r"$\Lambda$")
plt.ylim(0, 1000)

plt.tight_layout()
plt.subplots_adjust(hspace=0.5)
plt.show()
