
from obspy.signal import cross_correlation

#small box

bheight = 350*1000
bce = bheight*np.array([0,0,2,2,0])
bcn = bheight*np.array([0,1,1,0,0])
ang = np.deg2rad(-35)
rot = np.array([[np.cos(ang),-np.sin(ang)],[np.sin(ang), np.cos(ang)]])
bcer, bcnr = rot @ np.array([bce,bcn])

lre, lrn = proj.transform_point(-121.5,35.0, geo)

bcer += lre
bcnr += lrn

# even smaller box

# bheight = 150*1000
# bce = bheight*np.array([0,0,2,2,0])
# bcn = bheight*np.array([0,1,1,0,0])
# ang = np.deg2rad(-20)
# rot = np.array([[np.cos(ang),-np.sin(ang)],[np.sin(ang), np.cos(ang)]])
# bcer, bcnr = rot @ np.array([bce,bcn])

# lre, lrn = proj.transform_point(-119.,33.75, geo)

# bcer += lre
# bcnr += lrn


def utmtobox(e,n):
    e -= lre
    n -= lrn
    e, n = rot.T @ np.array([e,n])
    return (e, n)

def boxtoutm(e,n):
    e, n = rot @ np.array([e,n])
    return (e+lre, n+lrn)

for tr in data:
    e = tr.stats.easting
    n = tr.stats.northing
    bce, bcn = utmtobox(e,n)
    if (bce < 0) or (bce > bheight*2) or (bcn < 0) or (bcn > bheight):
        data.remove(tr)
    else:
        tr.stats.bce = bce
        tr.stats.bcn = bcn

# Quality Control
for tr in data:
    tr.stats.mlamp = np.mean(np.log(np.abs(hilbert(tr.data))))

reflamp = np.median(np.array([tr.stats.mlamp for tr in data]))
iqrlamp = (np.quantile(np.array([tr.stats.mlamp for tr in data]),0.75)-
           np.quantile(np.array([tr.stats.mlamp for tr in data]),0.25))

for tr1 in data:
    cc = np.array([cross_correlation.xcorr_max(cross_correlation.correlate(tr1, tr2, 500))[1] for tr2 in data]) #gets out cross correlation value
    p = np.linspace(0,100,101)
    qdf = np.percentile(np.abs(cc),p)
    tr1.stats.qdf = qdf 
    tr1.stats.goodflag = ((tr1.stats.qdf[51] > 0.60) and (np.abs(tr1.stats.mlamp-reflamp) < 3*iqrlamp))# this says that (100-50) = 50% of abs cc have to be better than 0.60 to be good, and that the rms log amplitude must be within 3 std dev of the average
        
data_reject = obspy.Stream()
for tr in data:
    if not tr.stats.goodflag:
        data_reject += tr
        data.remove(tr)