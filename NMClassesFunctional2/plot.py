import numpy as np
import math as m
import matplotlib.pyplot as plt

def f(x,V):
	match V[0]:
		case 0:
			return x - 1
		case 1:
			return x+m.exp(-x)
		case 2:
			return 1
		case 3:
			return (25 + 27 * m.cos(2*x))/(160 * m.pi)
		case 4:
			return (x+((6*x-2) * V[1] - V[1] * V[1])/( V[1] * V[1] - 3 * V[1] + 6))/6
		case 5:
			return m.cos(2*x)
		case 6:
			return 3 * x *(2 * V[1] - 3 * V[1] * x + 6) / (V[1]*V[1] - 18 * V[1] + 18)
		case 7:
			return 1 + 4 * x/9
		case 8:
			return x
		case 9:
			return 1
		case 10:
			return 1 + 2 * V[1] * m.cos(x) * m.cos(x) / (2 - V[1] * m.pi)
			


tit = [(r'$10.\/u(x) - \int_{0}^{1} u(t)dt = x - \frac{1}{2}, \/ u(x) = x - 1$'),
(r'$11.\/u(x) - \frac{1}{2} \int_{0}^{1} xe^{t}\/u(t)dt = e^{-x}, \/ u(x) = x + e^{-x}$'),
(r'$12.\/u(x) - \int_{0}^{0.5} \sin(xt)\/u(t)dt = 1 + \frac{1}{x}(\cos(\frac{x}{2})-1), \/  u(x) \equiv 1$'),
(r'$13.\/u(x) + \frac{1}{4\pi} \int_{0}^{2\pi}\frac{1}{\sin^{2}(\frac{x+t}{2}) + 0.25\cos^{2}(\frac{x+t}{2})}\/u(t)dt = \frac{1}{16\pi}(5+3\cos(2x)),\/ u(x)=\frac{1}{160\pi}(25+27\cos(2x))$'),
(r'$14.\/u(x) - \lambda \int_{0}^{1}(2x-t)\/u(t)dt=\frac{x}{6}, \/ u(x)=\frac{1}{6}(x + \frac{(6x-2)\lambda - \lambda^{2}}{\lambda^{2} - 3 \lambda + 6})$'),
(r'$15.\/u(x) - \int_{0}^{2\pi}\sin(x)\cos(x)\/u(t)dt = \cos(2x), \/ u(x) = \cos(2x)$'),
(r'$16.\/u(x) - \lambda \int_{0}^{1}(4xt-t^{2})\/u(t)dt = x, \/ u(x) = \frac{3x(2\lambda - 3 \lambda x + 6)}{\lambda^{2} - 18\lambda + 18}$'),
(r'$17.\/u(x) - \int_{0}^{1}xt^{2}\/u(t)dt = 1, \/ u(x) = 1 + \frac{4}{9}x$'),
(r'$18.\/u(x) - \frac{1}{2} \int_{0}^{1}xt\/u(t)dt = \frac{5}{6}x, \/ u(x) = x$'),
(r'$19.\/u(x) - \int_{-1}^{1}x^{2}e^{xt}\/u(t)dt = 1 - x(e^{x} - e^{-x}), \/ u(x) \equiv 1$'),
(r'$20.\/u(x) - \lambda \int_{0}^{\pi}\cos^{2}(x)\/u(t)dt = 1, \/ u(x) = 1 + \frac{2\lambda}{2 - \pi\lambda}\cos^{2}(x)$')]

X = np.genfromtxt("X.txt", delimiter=",")
U = np.genfromtxt("U.txt", delimiter=",")
V = np.genfromtxt("V.txt")

print(V[0])
F = [f(x,V) for x in X]


plt.rcParams["mathtext.fontset"] = "cm"
fig, ax = plt.subplots(figsize = (10,8))

##fig, ax = plt.subplots()

ax.plot(X,U, label = r'$u_{approx}(x)$', color = 'C1', marker='o', markersize = 1)
ax.plot(X,F, label = r'$u(x)$',color ='C0')
ax.set_title(r'Fredholm integral equation of the second kind solved by @SAristeev' + "\n" + tit[int(V[0])], fontsize = 15, pad = 15)
ax.set_ylabel(r"$u(x)$", fontsize=15)
ax.set_xlabel(r"$x$" + "\nIntegral = " + str(V[2]) + " Number of nodes = " + str(int(V[3])), fontsize=15)
ax.grid(True)

ax.legend(prop={'size': 15})
plt.show()