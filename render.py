import pylab
from scipy import array, zeros, average, tile, ones, sqrt, sum, pi, r_, c_, \
    dot, linspace, sin, cos, exp, log, rand, interp, floor, cumsum, arange
from scipy.linalg import eig, norm
from enthought.mayavi import mlab
import asym_rotor as ar

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

mau = 196.97

echarge = 1.602e-19

def center_of_mass(q, m):
    v = zeros(3)
    for i in range(3):
        v[i] = sum(q[:,i]*m[:])/sum(m[:])
    return v

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

au9acm = au9a - tile(center_of_mass(au9a, mau9a), (au9a.shape[0], 1))
au9bcm = au9b - tile(center_of_mass(au9b, mau9b), (au9b.shape[0], 1))

au9adipole = array([ -2.78378518e-02,   2.08261165e-01,   2.18359268e-08])
#au9bdipole = array([  6.66498076e-02,   1.93246771e-10,   6.42467022e-11])
au9bdipole = array([  0.1,   0.0,   0.0])


linrot = array([
        [ 0, 0, 0.5],
        [ 0, 0, -0.5]])
m_linrot = array([5 * mau, 5 * mau])
linrot_dipole = array([0, 0, 4])

octrot = array([
[0, 0, 0],
[-2, 0, 0],
[2, 0, 0],
[0, 2, 0],
[0, -2, 0],
[0, 0, 1],
[0, 0, -1]])
m_octrot = 10.0 * mau * ones(len(octrot))/len(octrot)
octrot_dipole = array([0, 0, 4.0])

# UNITS:
# coords in angstrom
# mass in amu
# time in ps
# charge in e = 1.6e-19 C
# electric field in amu * angstrom / ps**2.0 / e = 1.04e6 V/m
# energy in amu * angstrom**2.0 / ps**2.0 = 1.67e-23 J = 0.1 meV
# temperature in Kelvin kB = 0.82

kB_au = 0.82

# a nice choice for field 19000 kV / 2.6 mm = F = 7.02

def inertia_tensor(q, m):
    I = zeros((3,3))
    for (i,j) in [(a,b) for a in range(3) for b in range(3)]:
        if i == j:
            I[i,j] = sum(m[:] * (q[:,0]**2.0 + q[:,1]**2.0 + q[:,2]**2.0 - 
                                 q[:,i]*q[:,j]))
        else:
            I[i,j] = - sum(m[:]*q[:,i]*q[:,j])
    return I

    
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
def total_energy(qt, F, mu, mol_I):
    """ compute the total energy """
    t = qt[0]
    q = qt[1:5]
    w = qt[5:8]
  
    Q = q_to_rotmatrix(q)
    T = dot(mol_I, w**2.0)/2.0 # kinetic energy
    V = dot(mu, dot(Q.T, array([0.0, 0.0, F]))) # potential energy of dipole
    
    return (T - V) * kB_au

def angular_momentum(qt, F, mu, mol_I):
    t = qt[0]
    q = qt[1:5]
    w = qt[5:8]

    Q = q_to_rotmatrix(q)

    return dot(Q.T, mol_I*w)

def dipole_projection(qt, F, mu, mol_I):
    """ 
    compute projection of mu onto the F vector in the space fixed
    axes.  This is is proportional to the instantaneous force on the particle
    """
    t = qt[0]
    q = qt[1:5]
    w = qt[5:8]
    
    Q = q_to_rotmatrix(q)
    
    return dot(array([0.0, 0.0, 1.0]), dot(Q, mu))


def initial_cond(coords, mass, dipole, temp, F):
    cm_coords = coords - tile(center_of_mass(coords, mass), 
                              (coords.shape[0], 1))

    print "computing inertia tensor and principal axes of inertia"

    mol_I, mol_Ix = eig(inertia_tensor(cm_coords, mass))
    mol_I.sort()

    print "principal moments of inertia are: ", mol_I
    
    # compute the ratio of the dipole energy to the 
    # rotational energy

    print "x = (mu*F / kB*T_R) = ", norm(dipole) * F / kB_au / temp

    # random initial angular velocity vector
    # magnitude set so that 0.5 * I * w**2.0 = kT 
    w_mag = sqrt(2.0 * kB_au * temp / mol_I.mean())
    w0 = 2.0*rand(3) - 1.0
    w0 = w0 / norm(w0) * w_mag

    # random initial orientation / random unit quaternion
    q0 = 2.0*rand(4) - 1.0
    q0 = q0 / norm(q0)
    
    return q0, w0

def interpolate_quat(t,qt,steps):
    """
    interpolate a sequence of quaternions to constant time

    intended only to make the animation smooth
    this is possibly a numerically dodgy procedure
    real calculation will use the actual runge-kutta timesteps
    """
    ti = linspace(min(t), max(t), steps)

    # interpolate the steps

    qti = zeros((steps, 4))
    qti[:,0] = interp(ti, t, qt[:,0])
    qti[:,1] = interp(ti, t, qt[:,1])
    qti[:,2] = interp(ti, t, qt[:,2])
    qti[:,3] = interp(ti, t, qt[:,3])
    
    # normalize the quaternion
    for i in xrange(steps):
        qti[i,:] /= norm(qti[i,:])
    
    return ti, qti


def animate(coords, mass, dipole, temp, F, cycles, w0=None, max_tstep=100000,
		no_anim=True):
    """
    UNITS:
    coords    angstrom
    mass      amu
    t         ps
    dipole    e * angstrom
    temp      kelvin 
    F         amu * angstrom / ps**2.0 / e
    """

    cm_coords = coords - tile(center_of_mass(coords, mass), 
                              (coords.shape[0], 1))

    mol_I, mol_Ix = eig(inertia_tensor(cm_coords, mass))
    mol_I.sort()
       
    print "principal moments of inertia"
    print mol_I
    
    q0,w0 = initial_cond(coords, mass, dipole, temp, F)

    print "initial conditions"
    print " q0 = ", q0, " norm(q0) = ", norm(q0)
    print " w0 = ", w0, " norm(w0) = ", norm(w0)
    
    #q0 = array([0., 0., 1., 0.])
    #w0 = array([0.,1.,0.], dtype=float)

    print "time in picoseconds"

    # t0 = 2.0 * pi / norm(w0) * cycles[0]
    # t1 = t0 + 2.0 * pi / norm(w0) * cycles[1]
    # t2 = t1 + 2.0 * pi / norm(w0) * cycles[2]

    t0 = cycles[0]
    t1 = t0 + cycles[1]
    t2 = t1 + cycles[2]

    print "t0 = ", t0
    print "t1 = ", t1
    print "t2 = ", t2
    print "total time = ", t0 + t1 + t2

    print "integrating motion of rotor"

    qt = ar.asym_rotor(t0, t1, t2, max_tstep, 
                       r_[q0, w0], F, dipole, mol_I)

    print "qt.shape = ", qt.shape

    # now extract the orientation quaternions
    tsr = qt[:,0]
    qtr = qt[:,1:5]
    wtr = qt[:,5:8]

    # diagnostic plots

    # field in rk time steps
    ft = zeros(len(tsr))
    ft[tsr < t0] = 0.0;
    ft[(tsr > t0) * (tsr < t1)] = F * tsr / (t1 - t0)
    ft[tsr > t1] = F

    fp = pylab.figure(figsize=(3,4))
   	
    # total energy
    pylab.subplot(4,1,1)
    pylab.plot(tsr, array([total_energy(qt[i,:], ft[i], dipole, mol_I) 
                            for i in xrange(len(tsr))]))
    pylab.title('Total Energy')
    
    # angular momentum
    pylab.subplot(4,1,2)
    L = r_[[angular_momentum(qt[i,:], ft[i], dipole, mol_I) 
            for i in xrange(len(tsr))]]    
    
    print "L.shape = ", L.shape
    #pylab.plot(tsr, array([norm(L[i,:]) for i in xrange(len(tsr))]), 'r-')
    pylab.plot(tsr, array([L[i,0]for i in xrange(len(tsr))]), 'b-')
    pylab.plot(tsr, array([L[i,1]for i in xrange(len(tsr))]), 'g-')
    pylab.plot(tsr, array([L[i,2]for i in xrange(len(tsr))]), 'r-') 
    pylab.title('Angular Momentum')

	# dipole projection
    pylab.subplot(4,1,3)
    dp = array([dipole_projection(qt[i,:], ft[i], dipole, mol_I) 
                           for i in xrange(len(tsr))])
    dpavg = cumsum(dp)/(arange(len(dp))+1.0)
    pylab.plot(tsr, dpavg)
    

    pylab.title('Projected Dipole (polarization)') 
	
	# field strength
    pylab.subplot(4,1,4)
    pylab.plot(tsr, ft)
    pylab.title('Field Strength')

    pylab.subplots_adjust(hspace=0.95)

    pylab.show()

    if no_anim:
	return

    # frame rate is 10 ps to 1 s
    frames = floor(t2 / 10 * 30)
    print "total time is ", t2, " rendering ", frames, " frames"
    ts, qt = interpolate_quat(tsr, qtr, frames)

    # field
    ft = zeros(len(ts))
    ft[ts < t0] = 0.0;
    ft[(ts > t0) * (ts < t1)] = F * ts / (t1 - t0)
    ft[ts > t1] = F

    # render the molecule    
    f = mlab.figure(size=(500,500), bgcolor=(1,1,1), fgcolor=(0,0,0))
    # # draw the molecule
    molecule  = mlab.points3d(cm_coords[:,0], cm_coords[:,1], cm_coords[:,2],
                              scale_factor=2,
                              color=(1,1,0),
                              resolution=20)
    # and the dipole moment
    dip = mlab.quiver3d([0],[0],[0], 
                        [dipole[0]],[dipole[1]],[dipole[2]],
                        scale_factor=1,
                        line_width=2.0,
                        color=(0,1,0))
    
    mlab.view(distance=45.0)

    #timelab = mlab.text(0.1, 0.9, 'Time: 0 ps', width=0.4)
    #fieldlab = mlab.text(0.5, 0.9, 'Field: 0 kV/cm', width=0.4) 
    
    ms = molecule.mlab_source
    dms = dip.mlab_source
    for t, fi, q in zip(ts, ft, qt):
        #timelab.text = "Time: %d ps" % t
        #fieldlab.text = "Field: %f" % fi
        M = q_to_rotmatrix(q)
        cp = dot(cm_coords, M)
        dp = dot(dipole, M)*12.0
        ms.set(x = cp[:,0], y=cp[:,1], z=cp[:,2])
        dms.set(u = [dp[0]], v = [dp[1]], w = [dp[2]])


if __name__ == "__main__":
    # u = array([0.,1.,1.])
    # u = u / sqrt(sum(u**2.0))
    # const_rotation(au9b, ones(au9b.shape[0]),u)

    
    F = 0.0
    temp = 300.0
    cycles = (0,0,1000)
    #animate(au9b, mau9b, au9bdipole, temp, F, cycles)
    animate(octrot, m_octrot, octrot_dipole, temp, F, cycles, no_anim=False)
