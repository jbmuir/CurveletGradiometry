import numpy as np
from scipy.sparse import coo_matrix, bmat, identity
from scipy.sparse.linalg import eigs
import numpy.polynomial.polynomial as poly

def lysmer(nn, hv, vsv, vpv, rhov, f, modn = 1):
    """ Calculate the phase velocity, group velocity and mode shape of a Rayleigh wave

    Inputs:
        nn   -- number of nodes
        hv   -- vector of layer heights
        vsv  -- vector of shear velocities
        vpv  -- vector of compressional velocities
        rhov -- vector of densities
        f    -- frequency
        modn -- mode number (currently only fundamental), default = 1

    Outputs: (given in a dictionary)
        k    -- wavenumber at this frequency
        vpk  -- phase velocity at this frequency
        vgk  -- group velocity at this frequency
        edv  -- vertical displacement eigenfunction
        edh  -- horizontal displacement eigenfunction

    Written by:
        A translation of lysmer.m (Matt Haney, Victor Tsai) by Jack Muir
        August 2015

    Reference:
        Lysmer J. (1970). Lumped Mass Method for Rayleigh Waves.
        Bulletin of the Seismological Society of America, 60(1), 89–104.
    """
    #Make mu, lambda, omega
    muv     = rhov * vsv * vsv #mu = G
    lambdav = rhov * vpv * vpv - 2. * muv
    omega   = 2. * np.pi * f

    k1r, k1c, k1v = [], [], []
    k2r, k2c, k2v = [], [], []
    k3r, k3c, k3v = [], [], []
    mr, mc, mv    = [], [], []

    def p(l, i):
        return list(map(lambda x: x + 2 * i, l))
    for i in range(nn):
        #Alternative variables from Lysmer
        h       = hv[i]
        rho     = rhov[i]
        alpha   = ((2. * muv[i]) + lambdav[i]) / 6.
        beta    = muv[i] / 6.
        theta   = (muv[i] + lambdav[i]) / 4. #theta = phi
        psi     = (muv[i] - lambdav[i]) / 4.

        #form the stiffness and mass matrices
        if i != (nn-1):
            k1r += p([0, 2, 1, 3, 0, 2, 1, 3], i)
            k1c += p([0, 2, 1, 3, 2, 0, 3, 1], i)
            k1v += [2. * alpha * h] * 2 + [2. * beta * h] * 2 + [alpha * h] * 2 + [beta * h] * 2
            k2r += p([0, 1, 0, 3, 1, 2, 2, 3], i)
            k2c += p([1, 0, 3, 0, 2, 1, 3, 2], i)
            k2v += [2. * psi] * 2 + [2. * theta] * 2 +  [-2. * theta] * 2 + [-2. * psi] * 2
            k3r += p([0, 2, 1, 3, 0, 2, 1, 3], i)
            k3c += p([0, 2, 1, 3, 2, 0, 3, 1], i)
            k3v += [6. * beta / h] * 2 + [6. * alpha / h] * 2 + [-6. * beta / h] * 2 + [-6. * alpha / h] * 2
            mr  += p([0, 1, 2, 3], i)
            mc  += p([0, 1, 2, 3], i)
            mv  += [h * rho / 2.] * 4
        else:
            k1r += p([0, 1], i)
            k1c += p([0, 1], i)
            k1v += [2. * alpha * h] + [2. * beta * h]
            k2r += p([0, 1], i)
            k2c += p([1, 0], i)
            k2v += [2. * psi] * 2
            k3r += p([0, 1], i)
            k3c += p([0, 1], i)
            k3v += [6. * beta / h] + [6. * alpha / h]
            mr  += p([0, 1], i)
            mc  += p([0, 1], i)
            mv  += [h * rho / 2.] * 2

    k1 = coo_matrix((k1v, (k1r, k1c)), shape=(2 * nn, 2 * nn)).tocsr()
    k2 = coo_matrix((k2v, (k2r, k2c)), shape=(2 * nn, 2 * nn)).tocsr()
    k3 = coo_matrix((k3v, (k3r, k3c)), shape=(2 * nn, 2 * nn)).tocsr()
    m  = coo_matrix((mv, (mr,mc)), shape=(2 * nn, 2 * nn)).tocsr()

    # find the rayleigh wave speed which would exist
    # if the model were a halfspace with the minimum model velocity
    # this is a lower bound which can be passed to eigs

    vsmin = np.min(vsv)
    vpmin = np.min(vpv)
    t3 = 1. / vsmin ** 6.
    t2 = -8. / vsmin ** 4.
    t1 = 24. / vsmin ** 2. - 16. / vpmin ** 2.
    t0 = -16. * (1. - (vsmin / vpmin) ** 2.)
    roots  = poly.polyroots([t0, t1, t2, t3])
    rspeed = np.sqrt(np.min(roots))
    #set up the eigenvalue problem
    lhs = bmat([[None, identity(2 * nn, format="csr")], [omega * omega * m - k3, k2]])
    rhs = bmat([[identity(2 * nn, format="csr"), None], [None, k1]])
    #solve the eigenvalue problem
    w, v = eigs(lhs, k=modn, M=rhs, sigma=omega/rspeed, maxiter=50000)
    #the wavenumber
    k = w[0]
    #the phase velocity
    vpk = omega / k
    #the displacement eigenfunction:
    edun = v[:2*nn,0]
    #normalise the eigenfunction:
    factor = 1. / np.dot(edun, m.diagonal() * edun)
    ed = edun * np.sqrt(factor)
    # the group velocity
    darray = 2. * k * k1 - k2
    vgk = np.dot(edun, darray.dot(edun)) * factor / (2. * omega)

    return {"k" : k, "vpk" : vpk, "vgk" : vgk, "edh" : ed[::2], "edv" : ed[1::2]}
