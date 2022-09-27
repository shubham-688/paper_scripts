# this script finds the value of basement depth below a station
#using the drillhole depth data
# awk '{if ($3<0.595  && $3!="NaN") print($1,$2,$3*.68*1000,$4) }' P_delay_SA_all.txt > depth_SA_1.txt
# awk '{if ($3>0.595  && $3!="NaN") print($1,$2,$3*1.66*1000,$4) }' P_delay_SA_all.txt > depth_SA_2.txt
# cat depth_SA_*.txt > depth_SA_RF.txt

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
import sys
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Ellipse
#plt.ion()
####
def get_sub(x):
    normal = "Psb"
    sub_s = "ᵖˢᵇ"
    res = x.maketrans(''.join(normal), ''.join(sub_s))
    return x.translate(res)
 ###
def second_smallest(numbers):
  if (len(numbers)<2):
    return
  if ((len(numbers)==2)  and (numbers[0] == numbers[1]) ):
    return
  dup_items = set()
  uniq_items = []
  for x in numbers:
    if x not in dup_items:
      uniq_items.append(x)
      dup_items.add(x)
  uniq_items.sort()
  return  uniq_items[1]
####
# #######
# file=open('/Users/Shubham/Research/geological_data_straya/GSSA_basement/cover_thickness_2020.txt','r')
# lines=file.readlines()
# n_line=len(lines)
# x=[]
# y=[]
# z=[]
# for i in range(1,n_line):
#     split_line=lines[i].split()
#     xyz_t=[]
#     x.append(float(split_line[0].rstrip()))
#     y.append(float(split_line[1].rstrip()))
#     z.append(float(split_line[2].rstrip()))
#
# file.close()
# # coverTh=open("/Users/shubham/Research/Lake_eyre_data/geological_ft/drillholes_depthtobasement_shp/cover_thickness.txt","r")

#DISTANCE FUNCTION
def distance(x1,y1,x2,y2):
    d=np.sqrt((x1-x2)**2+(y1-y2)**2)
    return d
#
def idw_rblock(xz,yz,r,p):
    x_block=[]
    y_block=[]
    z_block=[]
    xr_min=xz-r
    xr_max=xz+r
    yr_min=yz-r
    yr_max=yz+r
    for i in range(len(x)):
        # condition to test if a point is within the block
        if ((x[i]>=xr_min and x[i]<=xr_max) and (y[i]>=yr_min and y[i]<=yr_max)):
            x_block.append(x[i])
            y_block.append(y[i])
            z_block.append(z[i])

    #calculate weight based on distance and p value
    w_list=[]
    for j in range(len(x_block)):
        d=distance(xz,yz,x_block[j],y_block[j]) #distance function is created outside this function
        if d>0:
            w=1/(d**p)
            w_list.append(w)
            z0=0
        else:
            w_list.append(0) #if meet this condition, it means d<=0, weight is set to 0

    #check if there is 0 in weight list
    w_check=0 in w_list
    if w_check==True:
        idx=w_list.index(0) # find index for weight=0
        z_idw=z_block[idx] # set the value to the current sample value
    else:
        wt=np.transpose(w_list)
        z_idw=np.dot(z_block,wt)/sum(w_list) # idw calculation using dot product
    return z_idw
#
# station_file=open("P_delay_SA_all_noNan.txt","r")
# station_basement_tps="station_basement_tps_idw_.5deg_p2.txt"
# f = open(station_basement_tps,'w')
# for line1 in station_file:
#     lat_st=round(float(line1.split()[0]),3)
#     long_st=round(float(line1.split()[1]),3)
#     RF_tps=float(line1.split()[2])
#     name= line1.split()[3]
#
#     radius=round(.5/2**.5,3) # half of block size (so we are looking for points within 0.5 degree)
#     # half block size = search radius / sqrt(2)
#     # so 0.5 degree serves as a block size of ~ 0.7
#     power=2 #
#     #A +ve number that defines the weight of the distance in the interpolation process. The weight of each known point decreases as the distance from it to the interpolated cell increases.
#     #The higher the value of the Power, the faster the weight of the known point decreases. If the Power is less than 1 the appearance of the resulting surface will be sharper.
#     #If the value of the Power is greater than 1 the appearance of the surface will be smoother.
#     #Very large values of the Power will result of a surface with only few known points influencing the value of the interpolated cell. The most commonly used value is 2.
#     depth_idw=idw_rblock(long_st,lat_st,radius,power)
#     # depth_idw_list.append(z_idw)
#
#     f.write('{} {} {} {} {}\n'.format(lat_st,long_st,RF_tps,name,round(depth_idw,3)))

#### this function calculates the Vs mean for a1 a2 at a given depth..
### so for depth=500m, output is mean of Vs (0 to 500m)
def Vs_avg_d(a1,a2,depth):
    d1=np.arange(0,depth+1,50)

    Vs=a1+a2*d1**.5
    return Vs.mean()
####################
drill=[]
std=[]
tps=[]
# dist_all=[]
st_name=[]
station_basement_tps="station_basement_tps_idw_.5deg_p2.txt"
for line3 in open(station_basement_tps,'r'):
    if float(line3.split()[2]) > .005 and line3.split()[4]!= 'nan':
        tps.append(float(line3.split()[2]))
        drill.append(float(line3.split()[4]))
        st_name.append(line3.split()[3])
        # number.append(float(line3.split()[6]))

######## Sorting lists based on tps
idx=np.argsort(np.array(tps))
tps1 = np.array(tps)[idx]
drill1 = np.array(drill)[idx]
st_name1 = np.array(st_name)[idx]
# dist_all1 = np.array(dist_all)[idx]
####################### interpolation deg3
# f3 = np.poly1d( np.polyfit(tps1, drill1, 3, full=False))#,w=1/(std1+0.000001)))
# a=np.polyfit(tps1, drill1, 3, full=True)#,w=1/(std1+0.000001))
# RMSE=float(np.sqrt(a[1]/len(tps1)))
# print('RMSE deg 3=',RMSE)
# ####################### interpolation expo
# fexp = np.poly1d( np.polyfit(tps1, np.log(drill1), 1, full=False))#,w=1/np.log(std1)))
# a=np.polyfit(tps1, np.log(drill1), 1, full=True)#,w=1/np.log(std1))
# RMSE=float(np.sqrt(a[1]/len(tps1)))
# print('RMSE exp=',RMSE)
# ####### divinding data in two sections
# drill_l=[]
# tps_l=[]
# dist_l=[]
# drill_h=[]
# tps_h=[]
# dist_h=[]
# for i in range(len(tps1)):
#     if tps1[i]> .595 and drill1[i]>500:
#         tps_h.append(tps1[i])
#         drill_h.append(drill1[i])
#         dist_h.append(dist_all1[i])
#
#     if tps1[i]< .61 and tps1[i]> .195 and drill1[i]<500:
#         tps_l.append(tps1[i])
#         drill_l.append(drill1[i])
#         dist_l.append(dist_all1[i])
# ######## Sorting lists based on tps
# idx=np.argsort(np.array(tps_l))
# tps_l1 = np.array(tps_l)[idx]
# drill_l1 = np.array(drill_l)[idx]
# dist_l1 = np.array(dist_l)[idx]
# ######## Sorting lists based on tps
# idx=np.argsort(np.array(tps_h))
# tps_h1 = np.array(tps_h)[idx]
# drill_h1 = np.array(drill_h)[idx]
# dist_h1 = np.array(dist_h)[idx]
# ##########
bins = np.linspace(0.0, 1.34, 16) # 16 for .088; 21 for 0.066; 13 for .1
digitized = np.digitize(tps1, bins)
bin_means_tps = np.array([tps1[digitized == i].mean() for i in range(1, len(bins))])
bin_means_tps = bin_means_tps[~np.isnan(bin_means_tps)]
bin_median_tps = np.array([np.median(tps1[digitized == i]) for i in range(1, len(bins))])
bin_median_tps = bin_median_tps[~np.isnan(bin_median_tps)]
bin_means_depth = np.array([drill1[digitized == i].mean() for i in range(1, len(bins))])
bin_means_depth= bin_means_depth[~np.isnan(bin_means_depth)]
bin_median_depth = np.array([np.median(drill1[digitized == i]) for i in range(1, len(bins))])
bin_median_depth = bin_median_depth[~np.isnan(bin_median_depth)]
bin_std_depth= np.array([drill1[digitized == i].std() for i in range(1, len(bins))])
bin_std_depth = bin_std_depth[~np.isnan(bin_std_depth)]
second_small=second_smallest(bin_std_depth)
# this step is done to make 0 std as non-zero by replacing them by next smallest std values.
for i in range(len(bin_std_depth)):
    if bin_std_depth[i] == 0:
        bin_std_depth[i]=second_small
#######
########## binning and finding mean and std for "divided data"
# bins = np.linspace(0.2, 0.6, 9) # 9 for .05; 11 for 0.04
# digitized = np.digitize(tps_l1, bins)
# bin_means_tps = [tps_l1[digitized == i].mean() for i in range(1, len(bins))]
# bin_means_depth = [drill_l1[digitized == i].mean() for i in range(1, len(bins))]
# bin_std_depth= [drill_l1[digitized == i].std() for i in range(1, len(bins))]
# ####### for t > 0.6
# bins_2 = np.linspace(0.6, 1.32, 9) # 10 for .08; 9 for .09; 7 for .12
# digitized_2 = np.digitize(tps_h1, bins_2)
# bin_means_tps_2 = np.array([tps_h1[digitized_2 == i].mean() for i in range(1, len(bins_2))])
# bin_means_tps_2=bin_means_tps_2[~np.isnan(bin_means_tps_2)]
# bin_means_depth_2 = np.array([drill_h1[digitized_2 == i].mean() for i in range(1, len(bins_2))])
# bin_means_depth_2= bin_means_depth_2[~np.isnan(bin_means_depth_2)]
# bin_std_depth_2= np.array([drill_h1[digitized_2 == i].std() for i in range(1, len(bins_2))])
# bin_std_depth_2=bin_std_depth_2[~np.isnan(bin_std_depth_2)]
##########
########## fitting two linear curves with binning bin_median_tps,bin_median_depth,yerr=bin_std_depth, dividing at tps=.58 bin

f_lin_l = np.poly1d( np.polyfit(bin_median_tps[:7], bin_median_depth[:7], 1, full=False))#,w=max(bin_std_depth[:7])/bin_std_depth[:7]))
a_l=np.polyfit(bin_median_tps[:7], bin_median_depth[:7], 1, full=True)#,w=max(bin_std_depth[:7])/bin_std_depth[:7])
RMSE_l=float(np.sqrt(a_l[1]/len(bin_median_tps[:7])))
print('RMSE deg 1 low=',RMSE_l)
#######
f_lin_h = np.poly1d( np.polyfit(bin_median_tps[6:], bin_median_depth[6:], 1, full=False,w=max(bin_std_depth[6:])/(bin_std_depth[6:])))
a_h=np.polyfit(bin_median_tps[6:], bin_median_depth[6:], 1, full=True,w=max(bin_std_depth[6:])/(bin_std_depth[6:]))
# a_h,c=np.polyfit(bin_means_tps_2, bin_means_depth_2, 1, cov=True,w=max(bin_std_depth_2)/(bin_std_depth_2))

RMSE_h=float(np.sqrt(a_h[1]/len(bin_median_tps[6:])))
print('RMSE deg 1 high=',RMSE_h)
#################################################### fitting two linear curves without binning, dividing at tps=.6
# f_lin_l = np.poly1d( np.polyfit(tps1[:159], drill1[:159], 1, full=False))#,w=max(dist_l1)/dist_l1))
# a_l=np.polyfit(tps1[:159], drill1[:159], 1, full=True)#,w=max(dist_l1)/dist_l1)
# RMSE_l=float(np.sqrt(a_l[1]/len(tps1[:159])))
# print('RMSE deg 1 low=',RMSE_l)
# #######
# f_lin_h = np.poly1d( np.polyfit(tps1[158:], drill1[158:], 1, full=False))#,w=max(dist_h1)/(dist_h1)))
# a_h=np.polyfit(tps1[158:], drill1[158:], 1, full=True)#,w=max(dist_h1)/(dist_h1))
# # a_h,c=np.polyfit(bin_means_tps_2, bin_means_depth_2, 1, cov=True,w=max(bin_std_depth_2)/(bin_std_depth_2))
#
# RMSE_h=float(np.sqrt(a_h[1]/len(tps1[158:])))
# print('RMSE deg 1 high=',RMSE_h)
####################################################
####################### interpolation deg2
# f2 = np.poly1d( np.polyfit(bin_median_tps, bin_median_depth, 2, full=False))#,w=max(bin_std_depth)/bin_std_depth))
# a=np.polyfit(bin_median_tps, bin_median_depth, 2, full=True)#,w=max(bin_std_depth)/bin_std_depth)
# RMSE=float(np.sqrt(a[1]/len(bin_median_tps)))
# print('RMSE deg 2=',RMSE)
####################################################

#### the following is to find the exp curve with the least RMSE
def rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())

# k1=np.linspace(100,1000,19)
# k1=np.linspace(200,600,41)
# # k2= np.linspace(1, 2, 21)
# k2= np.linspace(1.1, 1.9, 81)
# rmse_min=200
# rmse_min_all=[]
# x=bin_median_tps
# for i in range(0,len(k1)):
#     for j in range(0,len(k2)):
#
#         y=k1[i]*x*np.exp(x*k2[j]) #
#         rmse_val=rmse(y,bin_median_depth)
#         rmse_min_all.append(rmse_val)
#         if rmse_val < rmse_min:
#             rmse_min=rmse_val
#             i_min=i
#             j_min=j
# print("minimum RMSE = {}, for k1={} and k2={}".format(rmse_min,k1[i_min],k2[j_min] ))
#
# ##
# #### foloowing find minimum rmse curve for D=qt^2
# q=np.linspace(1500,1540,41)
# rmse_min=200
# rmse_min_all=[]
# x=bin_median_tps
# for i in range (len(q)):
#     y=q[i]*x**2
#     rmse_val=rmse(y,bin_median_depth)
#     rmse_min_all.append(rmse_val)
#     if rmse_val < rmse_min:
#         rmse_min=rmse_val
#         i_min=i
#
# print("minimum RMSE = {}, for q={} ".format(rmse_min,q[i_min] ))
######
# tps error from the file tpsb_error_stack_freq
bin_diff_4hz=[0.0167,0.04,0.1167,0.04,0.09,0.0825,0.07,0.08,0.025,0.01,0.014,0.02,0.0225]
bin_diff_4hz=np.array(bin_diff_4hz)
#################
carie_st=['GW22','KOOTA','YUDNA','GW25','GW24','GW26','GW27','GW28','GW30','GW31',\
'SOC09','WIRRA','GA05','GA07','ARCOO','PARAK','TL15','WPARC','OOAKD','SOGAP','BILLA','MTEBA']
#############
fig, ax=plt.subplots(figsize=(7,7))
# plt.style.use('seaborn-paper')
plt.style.use('seaborn')
# plt.style.use('seaborn-pastel')
# plt.style.use('ggplot')

x = np.linspace(.01, 1.5, 100)
x1 = np.linspace(.01, 0.58, 60)
x2 = np.linspace(.58, 1.5, 60)

rmse_min=200
rmse_min_all=[]
bin_tps=bin_median_tps
for a2 in range(26,42,2):
    for a1 in range(50,120,5):
        d=bin_median_depth
        d1=np.linspace(0,3000,61) # for plotting purposes

        Vs_avg1=np.zeros(61)
        for i in range(61):
            Vs_avg1[i]=Vs_avg_d(a1,a2,d1[i])

        n=len(bin_median_depth)
        Vs_avg=np.zeros(n)
        for i in range(n):
            Vs_avg[i]=Vs_avg_d(a1,a2,d[i])


        Vs = Vs_avg/1000
        Vs1 = Vs_avg1/1000

        Vp=.9409+(2.0947*Vs)-(0.8206*Vs**2)+(0.2683*Vs**3)-(0.0251*Vs**4)
        Vp1=.9409+(2.0947*Vs1)-(0.8206*Vs1**2)+(0.2683*Vs1**3)-(0.0251*Vs1**4)

        T=d*(Vp-Vs)/(Vp*Vs)/1000
        T1=d1*(Vp1-Vs1)/(Vp1*Vs1)/1000

        # d_bin=bin_tps*(Vp*Vs)*1000/(Vp-Vs)
        # plt.plot(T1,d1,'black',linewidth=.2,alpha=.45,zorder=0)

        rmse_val=rmse(T,bin_median_tps)
        rmse_min_all.append(rmse_val)
        if rmse_val < rmse_min:
            rmse_min=rmse_val
            a1_min=a1
            a2_min=a2

print("minimum RMSE = {}, for a1={} and a2={}".format(rmse_min,a1_min,a2_min ))
d=np.linspace(0,3000,61)

# Vs=a1_min+a2_min*d**.5
Vs_avg=np.zeros(61)
for i in range(61):
    Vs_avg[i]=Vs_avg_d(a1_min,a2_min,d[i])
Vs = Vs_avg/1000
Vp=.9409+(2.0947*Vs)-(0.8206*Vs**2)+(0.2683*Vs**3)-(0.0251*Vs**4)

T=d*(Vp-Vs)/(Vp*Vs)/1000
#plt.plot(T,d,'black',linewidth=1.25,alpha=.8,zorder=4)
# y=350*x*np.exp(x*k) # k= 1.4
# y=250*x*np.exp(x*k) # k= 1.75
# y=k1[i_min]*x*np.exp(x*k2[j_min]) # k= 1.75
#
# y2=1960*x**2-540*x # for exp
# y2=1508*x**2 # for quadratic

for i in range(len(drill1)):
    if st_name1[i] in carie_st:
        ax.scatter(tps1[i],drill1[i],c='MediumVioletRed',alpha=.75,zorder=2)#MediumVioletRed
    else:
        ax.scatter(tps1[i],drill1[i],c='CadetBlue',alpha=.8,zorder=2)
#unc
ax.plot(x1, f_lin_l(x1), color='Maroon', lw = 1.5,ls='--',zorder=4)
# unc
ax.plot(x1, f_lin_h(x1), color='Maroon', lw = 1.5,ls='--',alpha=.3,zorder=4)
# unc
ax.plot(x2, f_lin_l(x2), color='Maroon', lw = 1.5,ls='--',alpha=.3,zorder=4)
ax.plot(x2, f_lin_h(x2), color='Maroon', lw = 1.5,ls='--',zorder=4)

# ax.scatter(tps1,drill1,c='CadetBlue',alpha=.75)
# ax.plot(x, y, color='darkred', lw = 1.5,ls='--')
# ax.plot(x, f2(x), color='Maroon', lw = 1.5,ls='--') # for quad
# ax.plot(x, y2, color='Maroon', lw = 1.5,ls='--') # for quad but direct

# ax.errorbar(bin_means_tps,bin_means_depth,yerr=bin_std_depth,fmt='s',c='darkred',markersize=8,alpha=.5)
#unc
ax.errorbar(bin_median_tps,bin_median_depth,yerr=bin_std_depth,fmt='s',c='MidnightBlue',markersize=6,alpha=.75,zorder=6)
ax.errorbar(bin_median_tps,bin_median_depth,xerr=bin_diff_4hz,fmt='s',c='MidnightBlue',markersize=6,alpha=.75,zorder=6)

ax.errorbar(bin_median_tps[6],bin_median_depth[6],yerr=bin_std_depth[6],fmt='s',c='gold',markersize=6,alpha=.86,zorder=6)
ax.errorbar(bin_median_tps[6],bin_median_depth[6],xerr=bin_diff_4hz[6],fmt='s',c='gold',markersize=6,alpha=.86,zorder=6)

ax.grid(alpha=0.7)
# ax.scatter(tps1,drill1,c=np.array(number),cmap='cividis')

ax.tick_params(which='major',labelsize=13)
# ax.axvline(x=.2,ls='-',lw=.9,c='black',alpha=.45)
# ax.axvline(x=.67,ls='--',lw=.9,c='black',alpha=.45)

# ax.axhline(y=500,ls='--',lw=.9,c='black',alpha=.35)
ax.set_xlabel('Arrival Time Ps$_{b}$ (s)',fontsize=15)
ax.set_ylabel('Borehole basement depth (m)',fontsize=15)
#unc
# fig.text(0.5, 0.90, ' ', bbox={'facecolor': 'gray', 'alpha': 0.6, 'pad': 50})
fig.text(0.5, 0.94, 'D = {}T{} : T{} < 0.58 s'.format(round(a_l[0][0],1),'Ps$_{b}$','Ps$_{b}$'),fontsize=12,bbox={'facecolor': 'gray', 'alpha': 0.1, 'pad': 23}, color='Maroon',ha='center', va='center')
fig.text(0.5, 0.905, 'D = {}T{} {}: T{} ≥ 0.58 s'.format(round(a_h[0][0],1),'Ps$_{b}$',round(a_h[0][1],1),'Ps$_{b}$'),fontsize=12, color='Maroon',ha='center', va='center')
fig.text(0.5, 0.975, 'Best-fitting straight lines',fontsize=11, style='italic',weight='bold',color='Maroon',ha='center', va='center')
#########
# ax.annotate('D = {}T{} e^({}T{})'.format(k1[i_min],'$_{Psb}$',k2[j_min],'$_{Psb}$'), (.15,2750),color='Darkred',fontsize=12)
# ax.annotate('RMSE = {} m'.format(round(rmse_min,2)), (.15,2550),color='Darkred',fontsize=11)
# ax.annotate('D = {} T{}\u00b2'.format(round(q[i_min],1),'$_{Psb}$'), (.15,2750),color='Darkred',fontsize=12)
# ax.annotate('RMSE = {} m'.format(round(rmse_min,2)), (.15,2550),color='Darkred',fontsize=11)
#unc
#ax.annotate('RMSE = {} m'.format(round(RMSE_h,2)), (.85,2550),color='Maroon',fontsize=11)
#ax.annotate('RMSE = {} m'.format(round(RMSE_l,2)), (.85,2750),color='Maroon',fontsize=11)

#########

ax.scatter(.046,2230,s=80,c='MediumVioletRed',alpha=.85)
ax.scatter(.046,2060,s=80,c='CadetBlue',alpha=.85)
ax.annotate('Proterozoic Basins', (.085,2190),color='black',weight='bold',fontsize=11,alpha=.7)
ax.annotate('Phanerozoic Basins', (.085,2030),color='DarkSlateGray',weight='bold',fontsize=11)

#ax.annotate('Best-fitting velocity profile', (.046,2650),color='black',weight='bold',style='italic',fontsize=11)
#ax.annotate('Vs = 115 + 32d^.5', (.046,2500),color='black',fontsize=12)
# ax.annotate('Predicted T (sec)', (.085,2030),color='Darkgreen',fontsize=13)
# plt.plot(T*3.55,d,'black',alpha=.5)

ax.set_xbound(lower=-.05,upper=1.4)
ax.set_ybound(lower=-50,upper=3000)

plt.savefig('RF_tps_BH_IDW_.5_gji.pdf',bbox_inches='tight', pad_inches=0.2)


# ax.annotate('y = x+150', (3,150),family="cursive",color='Teal',fontsize=20)


fig.show()
