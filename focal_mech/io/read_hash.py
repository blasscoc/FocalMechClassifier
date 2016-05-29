"""
Read the demo data supplied by the HASH demo codes.
"""
import datetime
from numpy import genfromtxt, int_, deg2rad, array

def read_hash_solutions(filename):
    data = genfromtxt(filename)
    return dict(zip(int_(data[:,0]),deg2rad(data[:,21:24])))

def read_demo(phase_data, sta_reverse, reverse=True):
    """
    Read the demo data in and grab the polarity and stations coords.
    """
    
    polarity_data, event_time_map = parse_phase_file(phase_data)
    # Some of the stations were known to be recording the oppsite
    # polarity during certain periods: 
    reverse_data = parse_reverse(sta_reverse)

    if reverse:
        reverse_polarity(reverse_data, polarity_data, event_time_map)

    demo_data = {}
    for event, val in polarity_data.items():

        _, data = val        
        
        pol = array([deg2rad(data[-1]), deg2rad(data[-2]), data[1]])
        pol = pol.reshape([3,-1]).T

        demo_data[event] = pol

    return demo_data

def _parse_phase_header(line_iter):
    """  
    1-10 5i2 origin time, year, month, day, hour, minute 
    11-14 f4.2 origin time, seconds
    15-16 i2 latitude, degrees
    17 a1 latitude, 'S'=south
    18-21 f4.2 latitude, minutes
    22-24 i3 longitude, degrees
    25 a1 longitude, 'E'=east
    26-29 f4.2 longitude, minutes
    30-34 f5.2 depth, km
    35-36 f2.1 magnitude
    81-88 2f4.2 horizontal and vertical uncertainty, km
    """

    line = line_iter.next()

    cursor = 0

    year = int(line[cursor:cursor+2]); cursor += 2
    month = int(line[cursor:cursor+2]); cursor += 2
    day = int(line[cursor:cursor+2]); cursor += 2
    hour = int(line[cursor:cursor+2]); cursor += 2
    minute = int(line[cursor:cursor+2]); cursor += 2
    second = float(line[cursor:cursor+2]); cursor += 4

    event_id = int(line[122:138])

    second /= 100
    microsecond = 1.0E+6 * (second - int(second))

    # y2k bug:
    if year > 50:
        year = year + 1900
    else:
        year = year + 2000

    # seconds since epoch:
    t = (datetime.datetime(year, month, day, hour, minute, int(second),
                           int(microsecond)) - 
         datetime.datetime(1970,1,1)).total_seconds()

    return event_id, t
    
def _parse_polarity_data(line_iter):
    """
    1-4 a4 station name
    7 a1 polarity:U,u,+,D,d,or-
    8 i1 quality: 0=good quality, 1=lower quality, etc 
    59-62 f4.1 source-station distance (km)
    66-68 i3 takeoff angle
    79-81 i3  azimuth
    83-85 i3 takeoff angle uncertainty 
    87-89 i3 azimuth angle uncertainty 
    """

    station = []
    polarity = []
    quality = []
    src_sta_dist = []
    toa = []
    azimuth = []

    line = line_iter.next()
    while(line[0] != " "):
        cursor = 0


# From HASH driver1.f:
#30    continue
#        read (12,35) sname(k),pickpol(k),p_qual(k),
#     &           qdist,ith,iaz,isthe,isazi
#35      format (a4,2x,a1,i1,50x,f4.1,i3,10x,i3,1x,i3,1x,i3)
#
        station.append(line[cursor:cursor+4]); cursor += 6
        
        pol = line[cursor:cursor+1]; cursor += 1
        if pol.lower() == 'd':
            polarity.append(-1)
        else:
            polarity.append(1)
        
        quality.append(int(line[cursor:cursor+1])); cursor += 51
        src_sta_dist.append(float(line[cursor:cursor+4])/10); cursor += 4

        toa.append(float(line[cursor:cursor+3])); cursor += 13
        azimuth.append(float(line[cursor:cursor+3])); cursor += 3

        line = line_iter.next()

    event_id = int(line)

    return event_id, [station, polarity, quality, src_sta_dist, toa, azimuth]


def parse_phase_file(filename):
    """
    files consist of header/data
    """

    with open(filename,'r') as fp:
        lines = fp.readlines()

    polarity_data = {}
    event_time_map = {}
    line_iter = iter(lines)
    try:
        while(line_iter):
            ev_id, t = _parse_phase_header(line_iter)
            event_time_map[ev_id] = t
            ev_id, data = _parse_polarity_data(line_iter)
            polarity_data[ev_id] = (t, data)
    except:
        pass

    return polarity_data, event_time_map

def parse_reverse(filename):
    data = genfromtxt(filename,dtype=[('station','S4'),
                                      ('start_t', 'S8'),
                                      ('end_t', 'S8')])

    reverse = {}
    for dat in data:

        st = dat[1]

        if st != '0':
            year = int(st[:4])
            month = int(st[4:6])
            day = int(st[6:8])
        else:
            year = 1900
            month = 1
            day = 1
            
        tstart = (datetime.datetime(year, month, day) -
                     datetime.datetime(1970, 1, 1)).total_seconds()

        en = dat[2]

        if en != '0':
            year = int(en[:4])
            month = int(en[4:6])
            day = int(en[6:8])
        else:
            year = 2100
            month = 1
            day = 1

        tend = (datetime.datetime(year, month, day) -
                    datetime.datetime(1970, 1, 1)).total_seconds()
        
        reverse[dat[0]] = (tstart, tend)

    return reverse

def reverse_polarity(reverse_data, polarity, event_time_map):

    for key, val in polarity.items():
        t = event_time_map[key]

        t0 = val[0]
        stations = val[1][0]
        pol = val[1][1]

        updated_polarity = []
        for sta, p in zip(stations,pol):
            
            if sta in reverse_data.keys():

                if(t > reverse_data[sta][0] and
                        t < reverse_data[sta][1]):
                    #reverse the polarity:                
                    updated_polarity.append(p * -1)
                else:
                    updated_polarity.append(p)
        
            else:
                updated_polarity.append(p)

        #update polarity
        polarity[key][1][1] = updated_polarity
