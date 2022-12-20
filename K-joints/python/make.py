from splipy import *
from splipy.io import G2
import numpy as np
from numpy import cos, sin, sqrt, pi, arccos
from splipy.utils import rotation_matrix

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
    return (X,u,v)


def get_matching_shell(L,D,l,d,t,phi,offset,D2=None,d2=None):
    X,u,v   = cylinderIntersection(D/2, d/2, 41, phi)
    X2,u2,v2 = cylinderIntersection(D/2, (d+t)/2, 41, phi)
    if D2 and d2:
        X4,u4,v4 = cylinderIntersection(D2/2, d/2, 41, phi)
        crv4 = curve_factory.cubic_curve(X4.T, t=u4, boundary=curve_factory.Boundary.PERIODIC)
        crv4 = crv4.rebuild(p=3, n=26)
        X5,u5,v5 = cylinderIntersection(D/2, d2/2, 41, phi)
        crv5 = curve_factory.cubic_curve(X5.T, t=u5, boundary=curve_factory.Boundary.PERIODIC)
        crv5 = crv4.rebuild(p=3, n=26)

    tmax = max(v2)
    tmin = min(v2)

    # main intersection curve
    crv1 = curve_factory.cubic_curve(X.T, t=u, boundary=curve_factory.Boundary.PERIODIC)

    # weld intersection curve
    crv2 = curve_factory.cubic_curve(X2.T, t=u2, boundary=curve_factory.Boundary.PERIODIC)

    # end piece of chord
    crv3 = curve_factory.circle(r=d/2, center=(l*sin(phi),0,l*cos(phi)), normal=(sin(phi),0,cos(phi)), xaxis=(0,1,0))


    if True:
        crv1 = crv1.rebuild(p=3, n=26)
        crv2 = crv2.rebuild(p=3, n=26)
        crv3 = crv3.rebuild(p=3, n=26)

    srf = surface_factory.cylinder(r=D/2, h=2*L, center=(-L,0,0), axis=(1,0,0)).swap()
    srf2 = surface_factory.edge_curves(crv3, crv1)
    with G2('non-matching-shells.g2') as myfile:
        myfile.write(crv1)
        myfile.write(srf)
        myfile.write(srf2+(offset,0,0))
        myfile.write(srf2.mirror(normal=(1,0,0)).swap()-(offset,0,0))


    crv1 += (offset,0,0)
    crv2 += (offset,0,0)
    crv3 += (offset,0,0)
    crvsI = crv1.split([crv1.start(0) + i*(crv1.end(0)-crv1.start(0))/4.0 for i in range(4)])
    crvsW = crv2.split([crv2.start(0) + i*(crv2.end(0)-crv2.start(0))/4.0 for i in range(4)])
    crvsO = crv3.split([crv3.start(0) + i*(crv3.end(0)-crv3.start(0))/4.0 for i in range(4)])
#    for j in range(4):
#        if j%2 == 1: crvsW[j].reverse()
#        if j in [1,2]: axis = [-1,0,0]
#        else:          axis = [ 1,0,0]
#        for i in range(1,len(crvsW[j])):
#            R = rotation_matrix((len(crvsW[j])-i-1)/len(crvsW[j])*2*pi/100,axis)
#            crvsW[j][i,:] = crvsW[j][i,:] @ R
#        if j%2 == 1: crvsW[j].reverse()
    if D2 and d2:
        crv4 += (offset,0,0)
        crv5 += (offset,0,0)
        crvsI2 = crv4.split([crv4.start(0) + i*(crv4.end(0)-crv4.start(0))/4.0 for i in range(4)])
        crvsI3 = crv5.split([crv5.start(0) + i*(crv5.end(0)-crv5.start(0))/4.0 for i in range(4)])
        srfs1 = [surface_factory.edge_curves(c1,c2) for c1,c2 in zip(crvsI3, crvsW)]
        srfs2 = [surface_factory.edge_curves(c1,c2) for c1,c2 in zip(crvsO, crvsI2)]
    else:
        srfs1 = [surface_factory.edge_curves(c1,c2) for c1,c2 in zip(crvsI, crvsW)]
        srfs2 = [surface_factory.edge_curves(c1,c2) for c1,c2 in zip(crvsO, crvsI)]
    for s in srfs2:
        s.set_order(3,3)
        s.refine(0,7)
    srfs3 = [srf.clone().mirror(normal=(1,0,0)).swap() for srf in srfs1]
    srfs4 = [srf.clone().mirror(normal=(1,0,0)).swap() for srf in srfs2]

    crvL1 = crvsW[0].clone().mirror(normal=(1,0,0))
    crvL2 = crvsW[1].clone().mirror(normal=(1,0,0))

    srf5  = surface_factory.edge_curves(crvsW[0], crvL1)
    srf6  = surface_factory.edge_curves(crvsW[1], crvL2)
    srf5.set_order(3,3)
    srf6.set_order(3,3)
    srf5.refine(0,3)
    srf6.refine(0,3)
    chordEndTopRight1 = curve_factory.circle_segment((tmin-tmax)/2, D/2, center=(L,0,0), normal=(-L,0,0), xaxis=(0,cos(tmin),sin(tmin))).rebuild(p=3,n=8)
    chordEndTopRight2 = curve_factory.circle_segment((tmin-tmax)/2, D/2, center=(L,0,0), normal=(-L,0,0), xaxis=(0,0,1)).rebuild(p=3,n=8)
    srf7  = surface_factory.edge_curves(crvsW[2], chordEndTopRight2)
    srf8  = surface_factory.edge_curves(crvsW[3], chordEndTopRight1)
    srf7.set_order(3,3)
    srf8.set_order(3,3)
    srf7.refine(0,15)
    srf8.refine(0,15)
    chordEndBottomRight = curve_factory.circle_segment(2*pi-tmax+tmin, D/2, center=(L,0,0), normal=(L,0,0), xaxis=(0,-cos(tmin),sin(tmin))).rebuild(p=3,n=29)
    botCrvs = [chordEndBottomRight - (L-x,0,0) for x in [L,crvsW[0][0,0], crvL1[0,0], -L]]
    srfs9 = [surface_factory.edge_curves(c1,c2) for c1,c2 in zip(botCrvs[:-1],botCrvs[1:])]
    for s in srfs9:
        s.set_order(3,3)
    srfs9[0].refine(0,15)
    srfs9[1].refine(0,3)
    srfs9[2].refine(0,15)
    srfs10 = [srf.clone().mirror(normal=(1,0,0)).swap() for srf in [srf7,srf8]]

    return srfs1 + srfs2 + srfs3 + srfs4 + [srf5,srf6,srf7,srf8] + srfs9 + srfs10

def get_stupid_connecting_piece(D,H,d,h,offset):
    X1,u1,v1   = cylinderIntersection((D+H/2)/2, (d+h/2)/2, 41, phi)
    X2,u2,v2   = cylinderIntersection((D+H/2)/2, (d-h/2)/2, 41, phi)
    X3,u3,v3   = cylinderIntersection((D-H/2)/2, (d+h/2)/2, 41, phi)
    X4,u4,v4   = cylinderIntersection((D-H/2)/2, (d-h/2)/2, 41, phi)
    crv1 = curve_factory.cubic_curve(X1.T, t=u1, boundary=curve_factory.Boundary.PERIODIC).rebuild(p=3,n=26) + (offset,0,0)
    crv2 = curve_factory.cubic_curve(X2.T, t=u2, boundary=curve_factory.Boundary.PERIODIC).rebuild(p=3,n=26) + (offset,0,0)
    crv3 = curve_factory.cubic_curve(X3.T, t=u3, boundary=curve_factory.Boundary.PERIODIC).rebuild(p=3,n=26) + (offset,0,0)
    crv4 = curve_factory.cubic_curve(X4.T, t=u4, boundary=curve_factory.Boundary.PERIODIC).rebuild(p=3,n=26) + (offset,0,0)
    srf1 = surface_factory.edge_curves(crv1,crv2)
    srf2 = surface_factory.edge_curves(crv3,crv4)
    vol1 = volume_factory.edge_surfaces(srf1,srf2)
    vol1.set_order(3,3,3)
    vols = vol1.split([vol1.start(0) + i*(vol1.end(0)-vol1.start(0))/4.0 for i in range(4)])
    vols2 = [v.clone().mirror(normal=(1,0,0)).swap() for v in vols]
    return vols + vols2


### dimensions of chord
L = 10 # length
D = 2  # diameter
H = .15# thickness

### dimensions of chord
l = sqrt(50) # length
d = 1.5      # diameter
h = .1       # thickness

### other parameters
t = 0.3       # mesh weld width
offset = 2.5  # distance between joints
phi = pi/5    # angle of attack

outer = get_matching_shell(L,D+H/2,l,d+h/2,t,phi, offset)
inner = get_matching_shell(L,D-H/2,l,d-h/2,t,phi, offset, D+H/2, d+h/2)
vols = [volume_factory.edge_surfaces(srf2,srf1) for srf1,srf2 in zip(inner,outer)]
for v in vols:
    v.set_order(3,3,3)
# vols += get_stupid_connecting_piece(D,H,d,h,offset)
with G2('volumetric.g2') as myfile:
    myfile.write(vols)

with G2('matching-shells.g2') as myfile:
    myfile.write(outer)

with G2('all-shells.g2') as myfile:
    myfile.write(outer)
    myfile.write(inner)
