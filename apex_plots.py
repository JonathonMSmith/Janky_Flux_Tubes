import pysat
import datetime as dt
import matplotlib.pyplot as plt

print('making inst')
ivm = pysat.Instrument(platform='cnofs', name='ivm')

start = dt.datetime(2014,3,1)
stop = dt.datetime(2014,3,2)
date_array = pysat.utils.time.create_date_range(start, stop)

output_resolution = (320, 240)
my_dpi = 192
figsize = tuple([x/my_dpi for x in output_resolution])
save_dir = '/Users/jonathonsmith/apex_animation/frames'

print('starting iteration')
for date in date_array:
    ivm.load(date=date)
    if ivm.data.epmty:
        continue
    print('creating figure')
    fig, ax = plt.subplots(1,1, figsize=figsize)
    fig.figsize(figsize)
    ax.plot(ivm.data.index, ivm.data.apex_altitude)
    ax.set_lim(0, 900)
    ax.set_xlim(0, 24)
    ax.grid(True)
    fig.tight_layout()
    filename = os.join(save_dir, str(date))
    plt.savefig(filename, dpi=my_dpi)
    plt.close()
