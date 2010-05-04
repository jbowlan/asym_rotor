import scipy
import pylab
import numpy as nm
import ctypes as ct

from scipy import exp, sin, cos, tan, log, array, zeros, ones, r_, c_, dot, \
    pi 

import au

# these are some units used to measure sizes field strengths etc..

# useful physical constants
amu = 1.67e-27
debye = 3.33e-30
eps0 = 8.85e-12
fpe0 = 4. * pi * eps0
kB = 1.3806e-23
eV = 1.602e-19
ang = 1e-10

muB = 9.2741e-24 # J/T
mu0 = 4. * pi * 1e-7

# these are python routines to work with rigid body animations
# we want to convert quaternions into orientation vectors etc..

_libar = nm.ctypeslib.load_library('libar.so', '.')
_libar.asym_rotor.argtypes = [
    ct.c_double, ct.c_double,
    nm.ctypeslib.ndpointer(dtype = nm.double),
    nm.ctypeslib.ndpointer(dtype = nm.double),
    ct.c_ulong,
    ct.c_double,
    ct.c_double, ct.c_double, ct.c_double,
    ct.c_double, ct.c_double, ct.c_double
    ]

_libar.asym_rotor.restype = ct.c_int

def asym_rotor(t0, tf, max_tstep, qt0, F, mu, I):
    qt0 = nm.asarray(qt0, dtype=nm.double)
    if len(qt0) != 7:
        print "qt0 must have 7 arguments"
        return 
    
    # this should be a max_tstep x 8 array
    qt = nm.zeros((max_tstep,8), dtype=nm.double)

    mu1, mu2, mu3 = mu
    I1, I2, I3 = I
    
    ts = _libar.asym_rotor(t0, tf, 
                           qt0, 
                           qt, max_tstep, 
                           F, 
                           mu1, mu2, mu3, 
                           I1, I2, I3)
    
    print "ts = ", ts

    if ts < max_tstep:
        return qt[:ts,:]
    else: 
        return qt


# in all qt is a state variable
# qt = array([t, q0, q1, q2, q3, w1, w2, w3])

def q_to_rotmatrix(qt):
    """ map a quaternion to a rotation matrix """

    t = qt[0]
    q = qt[1:5]
    w = qt[5:8]

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
    print "running asym_rotor test"
    F = 1.0
    mu = [0.0, 0.0, 1.0]
    I = [1.0, 1.0, 1.0]
    
    
    qt = asym_rotor(0, 1000.0, 100000, 
                    [0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 10.0], 
                    0.1, 
                    [0.0, 0.0, 1.0], 
                    [1.0, 2.0, 3.0])
    print "qt.shape = ", qt.shape
    pylab.figure()
    pylab.plot(qt[:,0],qt[:,1])
    pylab.plot(qt[:,0],qt[:,2])
    pylab.plot(qt[:,0],qt[:,3])
    pylab.plot(qt[:,0],qt[:,4])
    pylab.show()
    












