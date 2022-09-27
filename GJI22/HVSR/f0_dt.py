#simple sciprt to plot point data

import matplotlib.pyplot as plt
import sys
### lat long dt(rf) DT(Tss) station_name f0 f0_std

fig, ax=plt.subplots(figsize=(6,6))
plt.style.use('seaborn-whitegrid')
# plt.style.use('bmh')



for line in open('5g_dt_DT_HVSR.txt','r'):
    line=line.split()
    if line[5]!='nan':
        dt=float(line[2])
        f0=float(line[5])
        std=float(line[6])
        txt=line[4]
        plt.errorbar(dt,f0, yerr=std, marker='o',fmt='None',color='darkgray',alpha=.85)
        ax.scatter(dt,f0, marker='o',color='cadetblue',zorder=10,alpha=.79)#darkkhaki
        # ax.annotate(txt, (dt+.02, f0),fontsize=7)
# for i, txt in enumerate(st_name):

ax.set_ylim(0.00001,1.2)
ax.set_xlim(0.00001,1.4)

plt.xlabel('Arrival Time Ps$_{b}$ (s)',fontsize=12)
plt.ylabel('f$_{0}$ (Hz)',fontsize=12)
plt.savefig('f0_tpsb.pdf',bbox_inches='tight', pad_inches=0.15)
plt.show()

#######################
sys.exit()
fig, ax=plt.subplots(figsize=(5,5))
plt.style.use('seaborn-whitegrid')

for line in open('5g_dt_DT_HVSR.txt','r'):
    line=line.split()
    if line[3]!='nan':
        tss=float(line[3])
        f0=float(line[5])
        t0=1/(2*f0)
        txt=line[4]
        # plt.errorbar(dt,f0, yerr=std, marker='o',fmt='None',color='darkkhaki',alpha=.95)
        ax.scatter(tss,t0, marker='o',color='cadetblue',zorder=10,alpha=.79)#


plt.plot((.45,1.55),(.45,1.55),'--', linewidth=1, markersize=1)
ax.set_ylim(.45,1.55)
ax.set_xlim(.45,1.55)
plt.xlabel('Tss from RF (sec)',fontsize=12)
plt.ylabel('T0 from HVSR',fontsize=12)
plt.savefig('T0_tss.pdf',bbox_inches='tight', pad_inches=0.15)
plt.show()
###
fig, ax=plt.subplots(figsize=(6,6))
plt.style.use('seaborn-whitegrid')
for line in open('5g_dt_DT_HVSR.txt','r'):
    line=line.split()
    if line[5]!='nan':
        dt=float(line[2])
        f0=float(line[5])
        std=float(line[6])
        txt=line[4]
        for line_ in open('tps_bh_5g.txt','r'):
            line_=line_.split()
            if line_[3]==line[4]:
                bh=float(line_[4])
                # plt.errorbar(bh,f0, yerr=std, marker='o',fmt='None',color='darkgray',alpha=.85)
                ax.scatter(bh,f0, marker='o',color='darkseagreen',zorder=10,alpha=.99)#darkkhaki

        # ax.annotate(txt, (dt+.02, f0),fontsize=7)
# for i, txt in enumerate(st_name):

ax.set_ylim(0.1,1.1)
ax.set_xlim(-50,3000)

plt.xlabel('Borehole basement depth (m)',fontsize=12)
plt.ylabel('f$_{0}$ (Hz)',fontsize=12)
plt.savefig('f0_bh.pdf',bbox_inches='tight', pad_inches=0.15)
plt.show()

##################
