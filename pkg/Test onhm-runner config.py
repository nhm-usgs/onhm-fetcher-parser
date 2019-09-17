import datetime
from fponhm import FpoNHM
# Determine the dates for the data pull
#start_pull_date, end_pull_date = compute_pull_dates(restart_date)
# DANGER!
start_pull_date = datetime.date(2019,9,1)
end_pull_date = datetime.date(2019, 9, 2)
print('pull period start = ', start_pull_date, ' end = ', end_pull_date)

# Run the Fetcher/Parser to pull available data
# RMCD: Note sure this works, I imagine this would be the call to make over-riding
# the START_DATE and END_DATE ENV variables setup in the nhmusgs-ofp Dockerfile
# The code below is one method to run ofp through nhumusgs/docker-images commented out
sformat = "%Y-%m-%d"
str_start_pull_date = start_pull_date.strftime(sformat)
str_end_pull_date = end_pull_date.strftime(sformat)
# # START_DATE and END_DATE are ENV variables in nhmusgs-ofp Docker file we over-ride with -e option
# ofp_docker_cmd = f"docker run ofp -e START_DATE={str_start_pull_date} END_DATE={str_end_pull_date}"
# with open("ofp_log.log", "a") as output:
#     subprocess.call(ofp_docker_cmd, stdout=output, stderr=output)
# Code below call fetch-parser thorough; assumes hru*.shp are in INDIR; will output to OUTDIR
INDIR = '../Data/'
OUTDIR = '../Output/'
print('starting Script')
# numdays = 2
fp = FpoNHM()
print('instantiated')
extract_type = 'days'
numdays = 2
# initialize(self, iptpath, optpath, weights_file, type=None, days=None, start_date=None, end_date=None)
ready = fp.initialize(INDIR, OUTDIR, INDIR + '/weights.csv', type='date', start_date=start_pull_date,
                      end_date=end_pull_date)
if ready:
    print('initalized\n')
    print('running')
    fp.run_weights()
    print('finished running')
    fp.finalize()
    print('finalized')
else:
    if extract_type == 'days':
        print('Gridmet not updated continue with numdays -1')
        fp.setNumdays(numdays - 1)
        print('initalized\n')
        print('running')
        fp.run_weights()
        print('finished running')
        fp.finalize()
        print('finalized')
    else:
        print('error: extract did not return period specified')