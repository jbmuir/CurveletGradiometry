import obspy
from obspy.clients.fdsn import Client

data_dir = "Data/"

client = Client("SCEDC", timeout=6000)

if not os.path.isfile(data_dir+"StationInventory/cistations.xml"):
    inventory = client.get_stations(network="CI",
                                   station="*",
                                   location="*",
                                   channel="BHZ",
                                   starttime = obspy.UTCDateTime(2010,1,1,0,0,0),
                                   endtime = obspy.UTCDateTime(2019,10,1,0,0,0),
                                   level="response")
    inventory.write(data_dir+"StationInventory/cistations.xml", format="STATIONXML")

inventory = obspy.read_inventory(data_dir+"StationInventory/cistations.xml")

evt_data_file = data_dir + "CMTLists/vanuatucmt.txt"

evlos, evlas = np.loadtxt(evt_data_file, usecols=(0,1), unpack=True)
evnames = np.loadtxt(evt_data_file, usecols=(12,), unpack = True, dtype=np.dtype('U13'))
evtimes = [obspy.UTCDateTime(int(e[0:4]), int(e[4:6]), int(e[6:8]), int(e[8:10]), int(e[10:12])) for e in evnames]


#event data
evname = evnames[data_no]
evlo = evlos[data_no]
evla = evlas[data_no]
evtime = evtimes[data_no]

data = obspy.Stream()

if not os.path.isfile(data_dir + "/Seismograms/{0}.mseed".format(evnames[data_no])):
    for sta in inventory[0].stations:
        stacode = sta.code
        channels = ["BHZ","BHE","BHN"]
        for chan in channels:
            try:
                data += client.get_waveforms(network="CI", 
                                station=stacode, 
                                location="*", 
                                channel=chan, 
                                starttime = evtimes[data_no], 
                                endtime = evtimes[data_no]+3600)
                print("Station {0} {1} success".format(stacode, chan))
            except:
                print("Station {0} {1} error".format(stacode, chan))
    data.write(data_dir + "/Seismograms/{0}.mseed".format(evnames[data_no]))

data = obspy.core.read(data_dir + "Seismograms/{0}.mseed".format(evnames[data_no]))

if only_BHZ:
    data = data.select(channel="BHZ")

if only_template:
    data = data.select(station="USC")

data.merge(fill_value=0)

data.detrend()

data.remove_response(inventory=inventory, output="DISP")

data.filter('bandpass', freqmin=1/100, freqmax=1/10, zerophase=True)

data.decimate(8, no_filter=True)

data.decimate(5, no_filter=True)
 
data.trim(starttime=evtime+2200, endtime=evtime+3000)

for (i, tr) in enumerate(data):
    network, station, location, channel = tr.id.split('.')
    chainfo = inventory.select(network=network,station=station,location=location,channel=channel)[0][0][0]
    tr.stats.latitude = chainfo.latitude
    tr.stats.longitude = chainfo.longitude
    easting, northing = proj.transform_point(tr.stats.longitude, tr.stats.latitude, geo)
    tr.stats.easting = easting
    tr.stats.northing = northing
    distance, az, baz = obspy.geodetics.gps2dist_azimuth(evla, evlo,tr.stats.latitude, tr.stats.longitude)
    tr.stats.distance = distance
