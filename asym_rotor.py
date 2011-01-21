import scipy as sp
import pylab as pl
import numpy as nm
import ctypes as ct

from scipy import exp, sin, cos, tan, log, array, zeros, ones, r_, c_, dot, \
    pi, newaxis, rand, randn, sqrt, flatnonzero, arange, cumsum

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

# in our nice unit system hbar = 6.311 amu A^2 / ps
hbar_au = 6.311
kB_au = 0.82


# these are python routines to work with rigid body animations
# we want to convert quaternions into orientation vectors etc..

_libar = nm.ctypeslib.load_library('libar.dylib', '.')
_libar.asym_rotor.argtypes = [
    ct.c_double, ct.c_double, ct.c_double,
    nm.ctypeslib.ndpointer(dtype = nm.double),
    nm.ctypeslib.ndpointer(dtype = nm.double),
    ct.c_ulong,
    ct.c_double,
    ct.c_double, ct.c_double, ct.c_double,
    ct.c_double, ct.c_double, ct.c_double
    ]

_libar.asym_rotor.restype = ct.c_int

_libar.asym_rotor_muz.argtypes = [
    ct.c_double, ct.c_double, ct.c_double,
    nm.ctypeslib.ndpointer(dtype = nm.double),
    ct.POINTER(ct.c_double),
    ct.c_ulong,
    ct.c_double,
    ct.c_double, ct.c_double, ct.c_double,
    ct.c_double, ct.c_double, ct.c_double
    ]


_libar.asym_rotor_muz.restype = ct.c_int

_libar.test_byref.argtypes = [ct.POINTER(ct.c_double)]
_libar.test_byref.restype = ct.c_int

def test_byref():
    t = ct.c_double()
    _libar.test_byref(ct.byref(t))
    return t


def asym_rotor(t0, t1, t2, max_tstep, qt0, F, mu, I):
    qt0 = nm.asarray(qt0, dtype=nm.double)
    if len(qt0) != 7:
        print "qt0 must have 7 arguments"
        return 
    
    # this should be a max_tstep x 8 array
    qt = nm.zeros((max_tstep,8), dtype=nm.double)

    mu1, mu2, mu3 = mu
    I1, I2, I3 = I
    
    ts = _libar.asym_rotor(t0, t1, t2, 
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

def asym_rotor_muz(t0, t1, t2, max_tstep, qt0, F, mu, I):
    qt0 = nm.asarray(qt0, dtype=nm.double)
    if len(qt0) != 7:
        print "qt0 must have 7 arguments"
        return 
   
    muz = ct.c_double()

    mu1, mu2, mu3 = mu
    I1, I2, I3 = I
    
    ts = _libar.asym_rotor_muz(t0, t1,t2, qt0, ct.byref(muz), max_tstep, F, mu1, 
                           mu2, mu3, I1, I2, I3) 
    
    print "ts = ", ts

    return muz

def inertia_tensor(q, m):
    I = zeros((3,3))
    for (i,j) in [(a,b) for a in range(3) for b in range(3)]:
        if i == j:
            I[i,j] = sum(m[:] * (q[:,0]**2.0 + q[:,1]**2.0 + q[:,2]**2.0 - 
                                 q[:,i]*q[:,j]))
        else:
            I[i,j] = - sum(m[:]*q[:,i]*q[:,j])
    return I

# routines to generate an array of random vectors
# making use of the built-in uniform and normal generators in numpy

def rand_vect_in_sphere(dim, n):
    qs = zeros((dim, n))
    rep = array(range(n))
    while(True):
        if len(rep) == 0:
            return qs
        qs[:,rep] = 2.0*rand(dim, len(rep)) - 1.0
        rep = rep[flatnonzero(sqrt(sum(qs[:,rep]**2.0, 0)) > 1.0)]

# I've been told that you can just take a nd normal vector and normalize it
# I would really like to take a few minutes to try and prove this however

def rand_unit_vect(dim, n):
    uvs = randn(dim,n)
    uvs /= sqrt(sum(uvs**2.0,0))
    return uvs


    
#def q_to_rotmatrix(q):
#    """ map a quaternion to a rotation matrix """
#    return array([[dot(q**2.0, array([+1.0, +1.0, -1.0, -1.0])),
#                   2.0*(q[1]*q[2] + q[0]*q[3]),
#                   2.0*(q[1]*q[3] - q[0]*q[2])],
#                  [2.0*(q[1]*q[2] - q[0]*q[3]),
#                   dot(q**2.0,array([+1.0, -1.0, +1.0, -1.0])),
#                   2.0*(q[2]*q[3] + q[0]*q[1])],
#                  [2.0*(q[1]*q[3] + q[0]*q[2]),
#                   2.0*(q[2]*q[3] - q[0]*q[1]),
#                   dot(q**2.0, array([+1.0, -1.0, -1.0, +1.0]))]])

def q_to_rotmatrix(q):
    """ map a quaternion to a rotation matrix """
    Q = array(
        [[  dot(q**2.0, array([+1.0, +1.0, -1.0, -1.0])),
            2.0*(q[:,1]*q[:,2] + q[:,0]*q[:,3]),
            2.0*(q[:,1]*q[:,3] - q[:,0]*q[:,2])],
         [  2.0*(q[:,1]*q[:,2] - q[:,0]*q[:,3]),
            dot(q**2.0,array([+1.0, -1.0, +1.0, -1.0])),
            2.0*(q[:,2]*q[:,3] + q[:,0]*q[:,1])],
         [  2.0*(q[:,1]*q[:,3] + q[:,0]*q[:,2]),
            2.0*(q[:,2]*q[:,3] - q[:,0]*q[:,1]),
            dot(q**2.0, array([+1.0, -1.0, -1.0, +1.0]))]])
    return Q.swapaxes(0,2).swapaxes(1,2)

# this should not change with time
def total_energy(qt, F, mu, mol_I):
    """ compute the total energy """
    t = qt[0]
    q = qt[1:5]
    w = qt[5:8]
  
    Q = q_to_rotmatrix(q)
    T = dot(mol_I, w**2.0)/2.0 # kinetic energy
    V = dot(mu, dot(Q.T, array([0.0, 0.0, F]))) # potential energy of dipole
    
    return (T - V) * kB_au

def muz_bar(qt, F, mu, mol_I):
    t = qt[:,0]
    q = qt[:,1:5]
    w = qt[:,5:8]
    
    Q = q_to_rotmatrix(q)

    # transform the dipole to space fixed coords
    fsp = dot(Q, mu)

    n = arange(fsp.shape[0])+1.0
    fspav = cumsum(fsp,axis=0) / n[:,newaxis]
    return fspav[-1,2]  



def stats(qt, F, mu, mol_I):
    """
    follow the motion of a trajectory
    
    compute
   
    kinetic energy
    potential energy
    total energy
    projected dipole
    total angular momentum
    Jx
    Jy
    Jz

    """
    t = qt[:,0]
    q = qt[:,1:5]
    w = qt[:,5:8]

    print "w.shape = ", w.shape
    print "q.shape = ", q.shape
    print "mol_I.shape = ", mol_I.shape

    # compute the third column of the inverse orientation matrix
    invQz = array([
        2.0 * (q[:,1]*q[:,3] + q[:,0]*q[:,2]),
        2.0 * (q[:,2]*q[:,3] - q[:,0]*q[:,1]),
        q[:,0]**2.0 - q[:,1]**2.0 - q[:,2]**2.0 + q[:,3]**2.0]) 


    # kinetic energy
    T = dot(w**2.0, mol_I)/2.0
    print "T.shape", T.shape

    V = dot(mu, F*invQz)
    print "V.shape", V.shape

    U = T + V
    
    # projected dipole
    Q = q_to_rotmatrix(q)
    print "Q.shape = ", Q.shape

    # transform the dipole to space fixed coords
    fsp = dot(Q, mu)
    n = arange(fsp.shape[0])+1.0
    fspav = cumsum(fsp,axis=0) / n[:,newaxis]
    
    print "fsp.shape = ", fsp.shape

    J = w * mol_I[newaxis,:]
    #print "Jbf.shape = ", Jbf.shape

    #J = zeros((len(t), 3))
    #J[:,0] = nm.sum(Q[:,0,:] * Jbf, axis=1)
    #J[:,1] = nm.sum(Q[:,1,:] * Jbf, axis=1)
    #J[:,2] = nm.sum(Q[:,2,:] * Jbf, axis=1)

    print "J.shape = ", J.shape

    return t, T, V, U, J[:,0], J[:,1], J[:,2], fsp, fspav



def angular_momentum(qt, F, mu, mol_I):
    t = qt[0]
    q = qt[1:5]
    w = qt[5:8]

    Q = q_to_rotmatrix(q)

    return dot(Q.T, mol_I*w)


# in all qt is a state variable
# qt = array([t, q0, q1, q2, q3, w1, w2, w3])

if __name__ == "__main__":
    print "running asym_rotor test"
    # parameters of PABN molecule

    au9I = array([  6385.54243039,   8817.03576255,  15202.57819294])
    au9mu = array([ 0.285, 0, 0 ]) * 0.208

    F = 60.10
    #mu = array([1.24, 0.0, 0.214])
    mu = array([1.0/sqrt(3), 1.0/sqrt(3), 1.0/sqrt(3)])

    #rot_consts = array([0.230, 0.0405, 0.034]) # given in amu A^2 / ps^2
    rot_consts = array([0.1, 0.1, 0.1]) # given in amu A^2 / ps^2
    
    mol_I = hbar_au / 2 / rot_consts
    
    mol_I = au9I
    mu = au9mu
    
    print "mol_I = ", mol_I

    #q0 = array([0.0, 1.0, 0.0, 0.0])
    q0 = rand_unit_vect(4,1).ravel()
    w0 = 0.2*rand_vect_in_sphere(3,1).ravel()

    t0 = 200.0e3
    t1 = 300.0e3
    t2 = 1000.0e3

    qt = asym_rotor(t0, t1, t2, 1000000, r_[q0, w0], F, mu, mol_I) 
           
    print "qt.shape = ", qt.shape
   
    print "I = ", mol_I

    print "kinetic energy = (K)", 0.5 * dot(w0.T, mol_I*w0) / kB_au
    print "pF / kT = ", sqrt(sum(mu**2.0)) * F / 0.5 / dot(w0.T, mol_I*w0)  
    print "time at final time step is (ns) = ", qt[-1,0] / 1000.

    print "now computing asym_rotor_muz: "
    muzcomp = asym_rotor_muz(t0, t1, t2, 1000000, r_[q0, w0], F, mu, mol_I)
    print "muzcomp = ", muzcomp

    tsr = qt[:,0]

    # field in rk time steps
    ft =  F * (tsr - t0) / (t1 - t0)
    ft[flatnonzero(tsr < t0)] = 0.0;
    ft[flatnonzero(tsr > t1)] = F

    t, T, V, U, Jx, Jy, Jz, fsp, fspav = stats(qt, ft, mu, mol_I)
    #t, T, V, U = stats(qt, ft, mu, mol_I)
    
    pl.figure(figsize=(4,4))
    
    pl.subplot(311)
    pl.plot(t/1e3, T, 'g-', lw=0.5, alpha=0.5, label="T")
    pl.plot(t/1e3, V, 'b-', lw=0.5, alpha=0.5, label="V")
    pl.plot(t/1e3, U, 'r-', lw=1.5, label="U")
    
    pl.subplot(312)
    #pl.plot(t/1e3, sqrt(nm.sum(qt[:,1:5]**2.0,axis=1)), 'r-', label="|q|")
    pl.plot(t/1e3, ft, 'r-', label="|q|")
    
    pl.subplot(313)
    #pl.plot(t/1e3, Jx, 'r-', alpha=0.5, label="Jx")
    #pl.plot(t/1e3, Jy, 'g-', alpha=0.5, label="Jy")
    #pl.plot(t/1e3, Jz, 'b-', alpha=0.5, label="Jz")
   
    #pl.subplot(414)
    pl.plot(t/1e3, fsp[:,2], 'b-', alpha=0.5, label="mu_z")
    pl.plot(t/1e3, fspav[:,2], 'r-', label="<mu_z>")
    
    
    #pl.plot(t/1e3, fsp[:,1], 'g-', label="mu_y")
    #pl.plot(t/1e3, fsp[:,2], 'b-', label="mu_z")
    pl.legend()
    pl.show()












