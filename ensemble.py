import scipy as sp
import pylab as pl
import numpy as nm
import ctypes as ct

from scipy import exp, sin, cos, tan, log, array, zeros, ones, r_, c_, dot, \
    pi, rand, sqrt, where, flatnonzero, randn

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

# in our nice unit system hbar = 6.311 amu A^2 / ps
hbar_au = 6.311
kB_au = 0.82


# some nice rotors

# spherical rotor
sphere_rot_consts = array([0.1, 0.1, 0.1])
sphere_mol_I = hbar_au / 2 / sphere_rot_consts
sphere_mu = array([1.0/sqrt(3), 1.0/sqrt(3), 1.0/sqrt(3)])

# PABN (para-aminobenzonitrile) molecule
pabn_mu = array([1.24, 0.0, 0.214])
pabn_rot_consts = array([0.230, 0.0405, 0.034]) # given in amu A^2 / ps^2
pabn_mol_I = hbar_au / 2 / pabn_rot_consts

#Au9 planar cluster
au9_mu = array([0.0, 0.0, 0.2])
au9_rot_consts = array([0.230, 0.0405, 0.034]) # given in amu A^2 / ps^2
au9_mol_I = hbar_au / 2 / au9_rot_consts

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


# generate a sequence of random orientations and velocities using the metropolis 
# method
def generate_ensemble(n, temp, F, mu, mol_I):
    qs = zeros((n, 7))
    
    """
    the angular velocity is chosen by applying the metropolis method
    for the initial distribution we choose a vector uniformly distributed inside 
    of a sphere of radius w_cf = A * sqrt( temp * kb_au / mean(mol_I) )
    where A is a cutoff factor, chosen so that the probability of encountering a 
    vector longer than w_cf is negligable
    """

    kT = kb_au * temp
    w_cf = A * sqrt( temp * kb_au / mean(mol_I) )

    qs = zeros((7,n))
    rep = array(range(n))
    while(True):
        if len(rep) == 0:
            return qs
        qs[0:4,rep] = rand_unit_vect(4, len(rep))
        qs[4:7,rep] = w_cf * rand_vect_in_sphere(3, len(rep))
        rep = rep[flatnonzero( rand(len(rep)) > exp(-ar.total_energy(qs[:,rep], 
            F, mu, mol_I)/kT))]


def compute_ensemble(qs, t0, t1, t2, F, mu, mol_I, max_steps = 1000000):
    muz = zeros(qs.shape[0])
    for i in xrange(qs.shape[0]):
        qt = ar.asym_rotor(t0, t1, t2, max_steps, qs[i,:], F, mu, mol_I)
        muz[i] = ar.muz_bar(qt,F,mu,mol_I)
        print "computing (", i, "/", qs.shape[0], ") <muz> = ", muz[i]
    return muz

def compute_profile(w_cf, F, mu, mol_I):
    qs = c_[rand_unit_vect(4,100).T, w_cf*rand_vect_in_sphere(3,100).T]
    muz = compute_ensemble(qs, 10e3, 20e3, 100e3, F, mu, mol_I)
    return muz
