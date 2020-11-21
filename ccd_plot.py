#!python3.8
import sys
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

filename = sys.argv[1]
filt = "R"

def read_name(filename):
	cospar = "none"
	norad = "none"
	name = "none"
	dt = "none"

	with open(filename) as myfile:
		head = [next(myfile) for x in range(10)]
	# print(head)
	# print("---------------------------------"
	for line in head:
		line = line.split()
		# print (line)
		if line[1].strip() == "COSPAR":
			cospar = line[3]
		if line[1].strip() == "NORAD":
			norad = line[3]
		if line[1].strip() == "NAME":
			name = " ".join(line[3:])
		if line[1].strip() == "dt":
			dt = line[3]
	return cospar, norad, name, dt


cospar, norad, name, dt = read_name(filename)
# print(cospar, norad, name, dt)



flux, mR, Az, El = np.loadtxt(filename, unpack=True, skiprows=11, usecols=(6, 8, 10, 11))
date, time = np.loadtxt(filename, unpack=True, skiprows=11, usecols=(0, 1), dtype={'names': ('date', 'time'), 'formats': ('S10', 'S12')})

# print (flux[0], mR[0])

date_time = []
for i in range(0, len(date)):
    date_time.append(datetime.strptime(date[i].decode('UTF-8') + ' ' + time[i].decode('UTF-8') + "000", "%Y-%m-%d %H:%M:%S.%f"))

# print (date_time[2])

## fig im MAG
plt.rcParams['figure.figsize'] = [12, 6]
dm = max(mR) - min(mR)
dm = dm * 0.1
plt.axis([min(date_time), max(date_time), max(mR) + dm , min(mR) - dm])
plt.plot(date_time, mR, "xr-", linewidth=0.5, fillstyle="none", markersize=3)


d, t = str(date_time[0]).split(" ")
plt.title("Satellite Name:%s, NORAD:%s, COSPAR:%s \n Date=%s  UT=%s   dt=%s  Filter=%s" % (name, norad, cospar, d, t, str(dt), filt), pad=8, fontsize=12)
plt.ylabel('m_st')
plt.xlabel('UT')
ax = plt.gca()

# Azimuth axis----------------------------------
ax2 = ax.twiny()
ax2.set_xlim(ax.get_xlim())
numElems = 5
tt_idx = np.round(np.linspace(0, len(date_time) - 1, numElems)).astype(int)
Tt2 = np.array(date_time)
Az2 = np.array(Az)
El2 = np.array(El)

# Az2s = ["%3.2f" % Azt for Azt in Az2]
Az2s = []
Tt2s = []
for kk in range(0, len(Az2)):
    azt = Az2[kk]
    elt = El2[kk]
    t = Tt2[kk]
    Az2s.append("%3.1f; %3.1f"%(azt, elt))
    Tt2s.append(t.strftime("%H:%M:%S"))
Az2s = np.array(Az2s)
ax2.set_xticks(Tt2[tt_idx])  # new_tick_locations
ax2.set_xticklabels(Az2s[tt_idx], fontsize=8)
ax2.set_xlabel(r"Az;h [deg]", fontsize=8)
ax2.tick_params(axis='x', which='major', pad=0)

Tt2s = np.array(Tt2s)
ax.set_xticks(Tt2[tt_idx])  # new_tick_locations
ax.set_xticklabels(Tt2s[tt_idx], fontsize=10)

ax.xaxis.grid()
ax.yaxis.grid()
# ----------------------------------------------------

# plt.show()
Tf = Tt2s[0]
Tf = Tf.split(":")
TF = Tf[0]+Tf[1]

plt.savefig(norad + "_UT" + TF + ".png")
plt.clf()
####


## fig im IMPULS-------------------------------------------------------------------
plt.rcParams['figure.figsize'] = [12, 6]
dm = max(flux) - min(flux)
dm = dm * 0.1
# print(flux[0])
plt.axis([min(date_time), max(date_time), min(flux) - dm, max(flux) + dm])
plt.plot(date_time, flux, "xr-", linewidth=0.5, fillstyle="none", markersize=3)


d, t = str(date_time[0]).split(" ")
plt.title("Satellite Name:%s, NORAD:%s, COSPAR:%s \n Date=%s  UT=%s   dt=%s  Filter=%s" % (name, norad, cospar, d, t, str(dt), filt), pad=8, fontsize=12)
plt.ylabel('Flux')
plt.xlabel('UT')
ax = plt.gca()

# Azimuth axis----------------------------------
ax2 = ax.twiny()
ax2.set_xlim(ax.get_xlim())
numElems = 5
tt_idx = np.round(np.linspace(0, len(date_time) - 1, numElems)).astype(int)
Tt2 = np.array(date_time)
Az2 = np.array(Az)
El2 = np.array(El)

# Az2s = ["%3.2f" % Azt for Azt in Az2]
Az2s = []
Tt2s = []
# flux2s = []
for kk in range(0, len(Az2)):
    azt = Az2[kk]
    elt = El2[kk]
    f = flux[kk]
    t = Tt2[kk]
    Az2s.append("%3.1f; %3.1f" % (azt, elt))
    Tt2s.append(t.strftime("%H:%M:%S"))
    # flux2s.append("%12.6f" % f)
Az2s = np.array(Az2s)
ax2.set_xticks(Tt2[tt_idx])  # new_tick_locations
ax2.set_xticklabels(Az2s[tt_idx], fontsize=8)
ax2.set_xlabel(r"Az;h [deg]", fontsize=8)
ax2.tick_params(axis='x', which='major', pad=0)

Tt2s = np.array(Tt2s)
ax.set_xticks(Tt2[tt_idx])  # new_tick_locations
ax.set_xticklabels(Tt2s[tt_idx], fontsize=10)


# flux2s = np.array(flux2s)
# ff_idx = np.round(np.linspace(0, len(flux) - 1, numElems)).astype(int)
# ax.set_yticks(flux[ff_idx])  # new_tick_locations
# ax.set_yticklabels(flux2s[ff_idx], fontsize=10)

import matplotlib.ticker as ticker
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%8.3f'))

ax.xaxis.grid()
ax.yaxis.grid()
# ----------------------------------------------------

# plt.show()
plt.savefig(norad + "_UT" + TF + "_Imp.png")
plt.clf()
####
# plt.savefig('line_plot.jpg', dpi=300, quality=80, optimize=True, progressive=True)  # 50% less size
