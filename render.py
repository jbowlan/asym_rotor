import pylab
from scipy import array, zeros, average, tile, ones, sqrt, sum, pi, r_, c_
from scipy.linalg import eig
from enthought.mayavi import mlab
import asym_rotor as ar

mau = 196.97

echarge = 1.602e-19

au9a = array([ 
    [     11.607089869793844,      5.996087010719368,      5.999999913697459],
    [     14.260846713870986,      6.062329540201776,      5.999999966813623],
    [      9.009614544414923,      6.005546950328713,      5.999999909783206],
    [      6.384618242752045,      6.177624483897441,      5.999999905632539],
    [     10.318473820400607,      8.296093760875156,      5.999999948054048],
    [     13.070398240880717,      8.322801503305650,      5.999999930700847],
    [      7.608618209793721,      8.410051646886625,      5.999999986294037],
    [     11.708782452264883,     10.533131210787859,      6.000000002192541],
    [      9.007902060603714,     10.601021807983226,      5.999999962328102]])


mau9a = mau * ones(au9a.shape[0])

au9b = array([
[     11.776449309591214,      8.399999999968557,      6.000000000000737],
[      8.994421130567343,      8.399999999938002,      5.999999999999682],
[      6.364903712058243,      8.399999999944480,      5.999999999998753],
[     10.303612658356220,     10.699447020208124,      6.000000000003402],
[     12.881997115778805,     10.751493933326458,      6.000000000005502],
[      7.684091744601174,     10.722164638134736,      6.000000000003002],
[     10.303612658376425,      6.100552979730839,      5.999999999997355],
[     12.881997115680623,      6.048506066664706,      5.999999999999072],
[      7.684091744621163,      6.077835361777611,      5.999999999997295]])

mau9b = mau * ones(au9b.shape[0])

def center_of_mass(q, m):
    v = zeros(3)
    for i in range(3):
        v[i] = sum(q[:,i]*m[:])/sum(m[:])
    return v

au9acm = au9a - tile(center_of_mass(au9a, mau9a), (au9a.shape[0], 1))
au9bcm = au9b - tile(center_of_mass(au9b, mau9b), (au9b.shape[0], 1))

au9adipole = array([ -2.78378518e-02,   2.08261165e-01,   2.18359268e-08])
au9bdipole = array([  6.66498076e-02,   1.93246771e-10,   6.42467022e-11])

def inertia_tensor(q, m):
    I = zeros((3,3))
    for (i,j) in [(a,b) for a in range(3) for b in range(3)]:
        if i == j:
            I[i,j] = sum(m[:] * (q[:,0]**2.0 + q[:,1]**2.0 + q[:,2]**2.0 - 
                                 q[:,i]*q[:,j]))
        else:
            I[i,j] = - sum(m[:]*q[:,i]*q[:,j])
    return I

def animate(coords, mass, dipole, q0, w0, F, t):
    # animate the motion of a molecule with molecular coordinates
    # initial orientation q0
    # angular velocity w0
    # for time t
    
    print "center the molecular coordinates on the center of mass"
    
    cm_coords = coords - tile(center_of_mass(coords, mass), 
                              (coords.shape[0], 1))

    print "computing inertia tensor and principal axes of inertia"
    
    mol_I, mol_Ix = eig(inertia_tensor(cm_coords, mass))

    mol_I.sort()

    print "principal moments of inertia are: ", mol_I

    print "integrating motion of rotor"

    qt = ar.asym_rotor(0.0, t, 100000, r_[q0, w0], F, dipole, mol_I)

    print "qt.shape = ", qt.shape

    return



    # now compute the trajectory by integrating the rigid body equations
    # of motion

    # f = mlab.figure(size=(500,500))
    
    # # draw the molecule
    # molecule  = mlab.points3d(au9bcm[:,0], au9bcm[:,1], au9bcm[:,2],
    #                           scale_factor=2,
    #                           color=(1,1,0),
    #                           resolution=20)
    
    # # and the dipole moment
    # dip = mlab.quiver3d([0],[0],[0], 
    #                     [au9bdipole[0]],[au9bdipole[1]],[au9bdipole[2]],
    #                     line_width=3.0,
    #                     color=(0,1,0)) 
    

    return p
