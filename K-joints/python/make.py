from splipy import *
from splipy.io import G2
import numpy as np
from numpy import cos, sin, sqrt, pi, arccos

def cylinderIntersection(R, r, n, phi=0, theta=0, sweep=(pi/2,5/2*pi)):
    """   returns sample points of the intersection curve between two cylinders.

    One major cylinder with larger radius which is along the x-axis, while a
    smaller cylinder with an angle of attack coming throught the large one.
    phi=theta=0 corresponds to the small cylinder along the y-axis.

    :param R: radius of the large cylinder
    :type R: float
    :param r: radius of the small cylinder (r<R)
    :type r: float
    :param n: number of evaluation points
    :type n: int
    :param phi: angle of attack in the xz-plane (along the major pipe). Range [0,pi/2)
    :type phi: float
    :param theta: angle of attack in the yz-plane (across the major pipe). Range: [0,pi/2)
    :type theta: float
    :param sweep: [start stop] angle sweep as seen from the smallest cylinder. Note that
                  start=0 refers to the points (x,y) = (1,0), and contiuing clockwise default [pi/2, 5/2*pi]
    :type sweep: tuple of two float
    :return: 3xn matrix of coordinates for n evaluation points
    :rtype: numpy.ndarray
    """

    assert( 0 <= theta <= pi/2)
    assert( 0 <= phi   <= pi/2)
    assert( r  < R)

    u = np.linspace(sweep[0], sweep[1], n)
    x = cos(phi) * r * cos(u) - sin(phi) * (-sin(phi) * r * cos(u) + sqrt(sin(phi) ** 2 * r ** 2 * cos(u) ** 2 - r ** 2 + R ** 2 + r ** 2 * cos(u) ** 2 * cos(phi) ** 2)) / cos(phi)
    y = -sin(theta) * sin(phi) * r * cos(u) + cos(theta) * r * sin(u) - sin(theta) * (-sin(phi) * r * cos(u) + sqrt(sin(phi) ** 2 * r ** 2 * cos(u) ** 2 - r ** 2 + R ** 2 + r ** 2 * cos(u) ** 2 * cos(phi) ** 2))
    z = R * sqrt(0.1e1 - (sin(theta) * sin(phi) * r * cos(u) - cos(theta) * r * sin(u) + sin(theta) * (-sin(phi) * r * cos(u) + sqrt(sin(phi) ** 2 * r ** 2 * cos(u) ** 2 - r ** 2 + R ** 2 + r ** 2 * cos(u) ** 2 * cos(phi) ** 2))) ** 2 / R ** 2)

    v  = pi - arccos((sin(theta) * sin(phi) * r * cos(u) - cos(theta) * r * sin(u) + sin(theta) * (-sin(phi) * r * cos(u) + sqrt(sin(phi)** 2 * r** 2 * cos(u)** 2 - r** 2 + R** 2 + r** 2 * cos(u)** 2 * cos(phi)** 2))) / R)

    X = np.vstack((x, y, z))
    return (X,v)


### dimensions of chord
L = 10 # length
D = 2  # diameter
H = .15# thickness

### dimensions of chord
l = sqrt(50) # length
d = 1        # diameter
h = .1       # thickness

### angle of attack
phi=pi/4

X,v = cylinderIntersection(D/2, d/2, 50, phi)

# main intersection curve
crv1 = curve_factory.cubic_curve(X.T, curve_factory.Boundary.PERIODIC)

# end piece of chord
crv2 = curve_factory.circle(r=d/2, center=(l*sin(phi),0,l*cos(phi)), normal=(sin(phi),0,cos(phi)), xaxis=(0,1,0))

srf = surface_factory.cylinder(r=D/2, h=2*L, center=(-L,0,0), axis=(1,0,0)).swap()
srf2 = surface_factory.edge_curves(crv2, crv1)
with G2('shells.g2') as myfile:
    myfile.write(crv1)
    myfile.write(srf)
    myfile.write(srf2+(D,0,0))
    myfile.write(srf2.mirror(normal=(1,0,0)).swap()-(D,0,0))
