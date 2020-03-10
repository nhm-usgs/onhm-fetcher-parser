import argparse
from pathlib import Path
import xarray as xr
import dask

def get_dm_files(idir, year):
    fprcp = None
    ftmax = None
    ftmin = None

    vprcp = 'prcp'
    vtmax = 'tmax'
    vtmin = 'tmin'

    fprcp = idir / vprcp / f'daymet_v3_{vprcp}_{year}_na.nc4'
    ftmax = idir / vtmax / f'daymet_v3_{vtmax}_{year}_na.nc4'
    ftmin = idir / vtmin / f'daymet_v3_{vtmin}_{year}_na.nc4'

    return xr.open_mfdataset([fprcp, ftmax, ftmin], combine='by_coords')

def main():
    idir = None
    odir = None
    dmyear = None
    wght_file = None

    my_parser = argparse.ArgumentParser(prog='Daymet_etl',
                                        description='Convert daymet to prms cbh in netcdf format')

    my_parser.add_argument('-i', '--inputdir', type=str,
                           help='Input directory containing daymet files',
                           default=None, required=True)
    my_parser.add_argument('-o', '--outputdir', type=str,
                           help='Directory to output netcdf prms cbh files',
                           default=None, required=True)
    my_parser.add_argument('-y', '--year', type=int,
                           help='Year of daymet file to permorm etl',
                           default=None, required=True)
    my_parser.add_argument('-w', '--weightsfile', type=str,
                           help='path/weight.csv - path to weight file', metavar='weight_file',
                           default=None, required=True)

    args = my_parser.parse_args()

    if args.inputdir is not None:
        idir = Path(args.inputdir)
        if not idir.exists():
            print(f'Input directory: {idir} - does not exist - exiting')
            return
    if args.outputdir is not None:
        odir = Path(args.outputdir)
        if not odir.exists():
            print(f'Output directory: {odir} - does not exist - exiting')
            return
    if args.year is not None:
        dmyear = args.year
    if args.weightsfile is not None:
        wght_file = Path(args.weightsfile)
        if not wght_file.exists():
            print(f'Weights file: {wght_file} - does not exist -exiting')

    dm_data = get_dm_files(idir, dmyear)

    print(dm_data)
if __name__ == "__main__":
    main()