import numpy as np
import matlab.engine

def create_basis(nx, ny):
    eng = matlab.engine.start_matlab()
    eng.addpath(r'.', nargout = 0)
    gmat, scales = eng.Curvelet_Basis_No_Dir(1.*nx, 1.*ny, nargout=2)
    G = np.array(gmat).T
    s = np.array(scales)
    gs = G.shape
    return (G.reshape(gs[0], nx, ny), s) #rather bizarrely, MATLAB requires non-integers here...

def create_bases():
    G, s = create_basis(16, 16)
    np.save('Basis/G_16_16.npy', G)
    np.save("Basis/s_16_16.npy", s)
    G, s = create_basis(32, 32)
    np.save('Basis/G_32_32.npy', G)
    np.save("Basis/s_32_32.npy", s)
    G, s = create_basis(32, 64)
    np.save('Basis/G_32_64.npy', G)
    np.save("Basis/s_32_64.npy", s)
    G, s = create_basis(64, 64)
    np.save('Basis/G_64_64.npy', G)
    np.save("Basis/s_64_64.npy", s)
    G, s = create_basis(64, 128)
    np.save('Basis/G_64_128.npy', G)
    np.save("Basis/s_64_128.npy", s)