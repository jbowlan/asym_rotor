import pylab
from scipy import array, zeros, average, tile, ones, sqrt, sum, pi, r_, c_, \
    dot, linspace, sin, cos, exp, log
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

# given coords, mass, and dipole
# in angstrom, amu, debye
def select_scale(coords, mass, dipole):
    cm_coords = coords - tile(center_of_mass(coords, mass), 
                              (coords.shape[0], 1))

    print "computing inertia tensor and principal axes of inertia"

    mol_I, mol_Ix = eig(inertia_tensor(cm_coords, mass))
    mol_I.sort()

    print "principal moments of inertia are: ", mol_I
    
    
    
    

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

    # now extract the orientation quaternions
    
    ts = qt[:,0]
    qt = qt[:,1:5]
    
    





def q_to_rotmatrix(q):
    """ map a quaternion to a rotation matrix """
    return array([[dot(q**2.0, array([+1.0, +1.0, -1.0, -1.0])),
                   2.0*(q[1]*q[2] + q[0]*q[3]),
                   2.0*(q[1]*q[3] - q[0]*q[2])],
                  [2.0*(q[1]*q[2] - q[0]*q[3]),
                   dot(q**2.0,array([+1.0, -1.0, +1.0, -1.0])),
                   2.0*(q[2]*q[3] + q[0]*q[1])],
                  [2.0*(q[1]*q[3] + q[0]*q[2]),
                   2.0*(q[2]*q[3] - q[0]*q[1]),
                   dot(q**2.0, array([+1.0, -1.0, -1.0, +1.0]))]])

# this should not change with time
def total_energy(qt, param):
    """ compute the total energy """
    t = qt[0]
    q = qt[1:5]
    w = qt[5:8]
    F = param[0]
    mu = param[1:4]
    I = param[4:7]

    Q = q_to_rotmatrix(qt)
    T = dot(I, w**2.0)/2.0 # kinetic energy
    V = dot(mu, dot(Q.T, array([0.0, 0.0, F]))) # potential energy of dipole
    
    return T - V

def angular_momentum(qt, param):
    t = qt[0]
    q = qt[1:5]
    w = qt[5:8]
    F = param[0]
    mu = param[1:4]
    I = param[4:7]

    # this vector should be the same in the space and body fixed coordinates
    return dot(I, w)

def dipole_projection(qt, param):
    """ 
    compute projection of mu onto the F vector in the space fixed
    axes.  This is is proportional to the instantaneous force on the particle
    """
    t = qt[0]
    q = qt[1:5]
    w = qt[5:8]
    F = param[0]
    mu = param[1:4]
    I = param[4:7]
    
    Q = q_to_rotmatrix(qt)
    
    return dot(array([0.0, 0.0, F]), dot(Q, mu))


if __name__ == "__main__":
    # u = array([0.,1.,1.])
    # u = u / sqrt(sum(u**2.0))
    # const_rotation(au9b, ones(au9b.shape[0]),u)

    mass = ones(au9b.shape[0])*197.0
    q0 = array([0.0, 0.0, 0.0, 1.0])
    w0 = array([0.0, 1.0, 0.0])
    F = 1e6
    t = 100.0
    animate(au9b, mass, au9bdipole, q0, w0, F, t)
            
            
def const_rotation(coords, mass, u):
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

    theta = linspace(0, 10*2*pi, 500)
    qt = [r_[cos(th/2),sin(th/2)*u] for th in theta] 

    # render the molecule

    f = mlab.figure(size=(500,500))
    
    # # draw the molecule
    molecule  = mlab.points3d(cm_coords[:,0], cm_coords[:,1], cm_coords[:,2],
                              scale_factor=2,
                              color=(1,1,0),
                              resolution=20)
    
    # and the dipole moment
    dip = mlab.quiver3d([0],[0],[0], 
                        [au9bdipole[0]],[au9bdipole[1]],[au9bdipole[2]],
                        scale_factor=8,
                        line_width=2.0,
                        color=(0,1,0))
    
    # mlab.title('Polar Au_9 Cluster')
    # mlab.axes(extent=[0, 5, 0, 5, 0, 5])


    ms = molecule.mlab_source
    dms = dip.mlab_source
    for q in qt:
        M = q_to_rotmatrix(q)
        cp = dot(cm_coords, M)
        dp = dot(au9bdipole, M)*12.0
        ms.set(x = cp[:,0], y=cp[:,1], z=cp[:,2])
        dms.set(u = [dp[0]], v = [dp[1]], w = [dp[2]])



