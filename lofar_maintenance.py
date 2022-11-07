import urllib
import time
import numpy as np
def getStation(line): return line.split("ddrivetip(")[1].split("'")[1].replace('<br>',',').split(',')

def getStarttime(line): return line.split('hideddrivetip()')[1].split('</TD><TD>')[1]
def getEndtime(line): return line.split('hideddrivetip()')[1].split('</TD><TD>')[2]
def getAllMaintenance(lofar_url="http://old.astron.nl/lofar-schedule/schedule/schedule.php"):
    """ Parse website that publishes lofar observations and maintenance. Obtain station and start and time for maintenance of core and remote stations """
    page = urllib.urlopen(lofar_url)
    maintenancelist=[(getStation(l),getStarttime(l),getEndtime(l)) for l in page if ('RESERVATION' in l or 'MAINTENANCE' in l) and ('RS' in l or 'CS' in l)]
    return maintenancelist
def getMaintenance(starttime,endtime,lofar_url="http://www.astron.nl/lofar-schedule/schedule/schedule.php"):
    """ Get all stations that are in maintenance between the starttime and endtime of the trigger window. Date format is 2018-05-09 09:00:00 """
    allMaintenance=getAllMaintenance(lofar_url)
    assert starttime < endtime, "starttime not before endtime"
    unavailable_stations=[]
    starttime=time.strptime(starttime,'%Y-%m-%d %H:%M:%S')
    endtime=time.strptime(endtime,'%Y-%m-%d %H:%M:%S')
    for stations,start,end in allMaintenance:
        start=time.strptime(start,'%Y-%m-%d %H:%M:%S')
        end=time.strptime(end,'%Y-%m-%d %H:%M:%S')
        if start <=starttime <= end:
            unavailable_stations.extend(stations)
        if start <= endtime <= end:
            unavailable_stations.extend(stations)
        if starttime <= start <= endtime:
            unavailable_stations.extend(stations)
        if starttime <= end <= endtime:
            unavailable_stations.extend(stations)

    unavailable_stations.extend('CS001')
    unavailable_stations.extend('CS032')

    return np.unique(unavailable_stations)
