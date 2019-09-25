import urllib3
import xmltodict
from datetime import datetime, timedelta
import sys

def getxml(url):
    http = urllib3.PoolManager()

    response = http.request('GET', url)
    try:
        data = xmltodict.parse(response.data)
    except:
        print("Failed to parse XML from response (%s)" % traceback.format_exc())
    return data


# "http://thredds.northwestknowledge.net:8080/thredds/ncss/grid/agg_met_tmmx_1979_CurrentYear_CONUS.nc/dataset.xml"
# http://thredds.northwestknowledge.net:8080/thredds/ncss/grid/agg_met_tmmn_1979_CurrentYear_CONUS.nc/dataset.xml
# http://thredds.northwestknowledge.net:8080/thredds/ncss/grid/agg_met_pr_1979_CurrentYear_CONUS.nc/dataset.xml
# http://thredds.northwestknowledge.net:8080/thredds/ncss/grid/agg_met_rmin_1979_CurrentYear_CONUS.nc/dataset.xml
# http://thredds.northwestknowledge.net:8080/thredds/ncss/grid/agg_met_rmax_1979_CurrentYear_CONUS.nc/dataset.xml
# http://thredds.northwestknowledge.net:8080/thredds/ncss/grid/agg_met_vs_1979_CurrentYear_CONUS.nc/dataset.xml

serverURL = 'http://thredds.northwestknowledge.net:8080/thredds/ncss/grid'

data_packets = ['agg_met_tmmn_1979_CurrentYear_CONUS.nc', 'agg_met_pr_1979_CurrentYear_CONUS.nc',
                'agg_met_rmin_1979_CurrentYear_CONUS.nc', 'agg_met_rmax_1979_CurrentYear_CONUS.nc',
                'agg_met_vs_1979_CurrentYear_CONUS.nc']
urlsuffix = 'dataset.xml'

now = datetime.today().date()
yesterday = now - timedelta(days=1)

for data in data_packets:
    masterURL = serverURL + '/' + data + '/' + urlsuffix
    datadef = getxml(masterURL)['gridDataset']['TimeSpan']['end']
    gm_date = datetime.strptime(datadef[:10],'%Y-%m-%d').date()
    if gm_date != yesterday:
        print(f'Gridmet data {data} is not available:\n' +
              'process exiting')
        sys.exit(1)
