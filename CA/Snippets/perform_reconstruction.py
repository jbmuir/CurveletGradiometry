import pywt
import celer
from tqdm import tqdm

wvt = 'db4'

data_wvt = np.array([np.hstack(pywt.wavedec(tr.data, wvt)) for tr in data]).T
ncoefs, nstations = data_wvt.shape


def fitcv(Gp, x):
    model = celer.LassoCV(cv=5, n_alphas=50, max_epochs=500000, fit_intercept=False)
    model.fit(Gp,x)
    return model.coef_, model.alpha_


if actually_do_reconstruction:
    results = []
    norms = []
    alphas = []
    alpha = 1.0

    os.makedirs("Results/{0}".format(evname), exist_ok=True)
    with tqdm(total=ncoefs) as pbar1:
        for dc in data_wvt:
            norm = np.std(dc)
            norms += [norm]
            f, a = fitcv(Gp, dc/norm)
            results += [f]
            alphas += [alpha]      
            pbar1.update(1)
            
    np.save("Results/{0}/fit{1}.npy".format(evname, channel), results)
    np.save("Results/{0}/norm{1}.npy".format(evname, channel), norms)
    np.save("Results/{0}/alphas{1}.npy".format(evname, channel), alphas)
    
    
if process_to_timedomain:
    zres = np.load("Results/{0}/fit{1}.npy".format(evname, "BHZ"))
    znorm  = np.load("Results/{0}/norm{1}.npy".format(evname, "BHZ"))
    zalpha = np.load("Results/{0}/alphas{1}.npy".format(evname, "BHZ"))
    zres *= znorm[:,np.newaxis]
    wvt_tmp = pywt.wavedec(data[0].data, wvt)
    wvt_lens = [len(wc) for wc in wvt_tmp]
    wvt_reconstructions_evZ = (zres@Gevp.T).T
    d_reconstructions_evZ = np.real(np.array([reconstruction(w, wvt_lens, wvt) for w in wvt_reconstructions_evZ]))
    wvt_reconstructions_Z = (zres@Gp.T).T
    d_reconstructions_Z = np.real(np.array([reconstruction(w, wvt_lens, wvt) for w in wvt_reconstructions_Z]))
    np.save("Results/{0}/td_{1}.npy".format(evname, "BHZ"), d_reconstructions_evZ)                             

