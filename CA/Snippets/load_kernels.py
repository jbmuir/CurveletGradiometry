from scipy.io import loadmat

if actually_load_kernels:
    crv = loadmat("../Curvelet_Basis_Construction/G_64_128.mat")
    G_mat = crv["G_mat"]
    crvscales = crv["scales"]
    G_matr = np.reshape(G_mat, (64,128,G_mat.shape[1]))
    
    #alternative if you have curvelab for matlab installed and have created your own kernels
#     G_mat = np.load("../Curvelet_Basis_Construction/Basis/G_64_128.npy")
#     crvscales = np.load("../Curvelet_Basis_Construction/Basis/s_64_128.npy")

bce = np.array([tr.stats.bce for tr in data.select(channel="BHZ")])
bcn = np.array([tr.stats.bcn for tr in data.select(channel="BHZ")])  
crv_bce = np.linspace(0,2*bheight,128)
crv_bcn = np.linspace(0,bheight,64)


if actually_load_kernels:
    cvtscaler = (2.0**(cscale*crvscales))
    G = []
    Gl = []
    for i in range(G_mat.shape[-1]):
        frame = rbs(crv_bce,crv_bcn,G_matr[:,:,i].T)
        G += [np.array(frame.ev(bce,bcn))]

    G = np.array(G).T
    Gn = np.std(G)
    G = G/Gn
    Gp = G/cvtscaler 
    Gevp = G_mat/cvtscaler  / Gn

crv_bce_mat, crv_bcn_mat = np.meshgrid(crv_bce, crv_bcn)
evbeastings = crv_bce_mat.flatten()
evbnorthings = crv_bcn_mat.flatten()
eveastings, evnorthings = boxtoutm(evbeastings,evbnorthings)
staeastings = np.array([tr.stats.easting for tr in data.select(channel="BHZ")])
stalats = np.array([tr.stats.latitude for tr in data.select(channel="BHZ")])
stanorthings = np.array([tr.stats.northing for tr in data.select(channel="BHZ")])
stalons = np.array([tr.stats.longitude for tr in data.select(channel="BHZ")])
evlon,evlat,_= ccrs.Geodetic().transform_points(proj,eveastings,evnorthings).T