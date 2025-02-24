#!/bin/csh
#ogr2ogr -f "GMT" aus5wgd_l.gmt  aus5wgd_l.shp -t_srs WGS84
set PS=map_oz.ps
#set Rmap=-Rg
#set Rmap=-R-135/40/-90/90
set Rmap=-R130/142/-32.5/-25 # Lake_eyre+ Flinders
set Rmap=-R117/152/-38/-9 # Lake_eyre+ Flinders

set Rmap=-R106.5/-43/152.5/-8r

#set J=-JM8/22
set J=-JM12/11 #sa
set J=-JS130/-30/4i #oz
#set J=-Jx0.023i/0.029i

# set topogrd = /Users/shubham/Research/etopo/etopo_sa_15s.grd
set topogrd = ~/Research/etopo/ETOPO1_Ice_g_gmt5.grd
# set topoxyz = ~/Research/etopo/etopo1.xyz.
set GAB = ~/Research/geological_data_straya/artesian_basin/81672/GAB_Hydrological_Boundary.gmt

set oleron = ~/Documents/ScientificColourMaps5/oleron/oleron.cpt
set tpushuf = ~/Documents/cpt-city/tp/tpushuf
set afrikakarte = ~/Documents/cpt-city/wkp/lilleskut/afrikakarte

gmt6 gmtset MAP_FRAME_TYPE plain MAP_FRAME_WIDTH 4p FONT_ANNOT_PRIMARY 9.5p MAP_TICK_LENGTH_PRIMARY 7p MAP_FRAME_PEN 2.2p
gmt6 gmtset COLOR_BACKGROUND lightsteelblue MAP_ANNOT_OBLIQUE 6 # MAP_ANNOT_OBLIQUE lat_horizontal

gmt6 psbasemap -BNEws -Bxa10f5 -Bya10f5 $J $Rmap -K >! $PS
gmt6 makecpt -C$afrikakarte -M -T-1620/890/20 > afrikakarte.cpt #  for transparency

###########
gmt6 grdimage $topogrd $J $Rmap -Bx -By -FRgray -Cafrikakarte.cpt -I+nt.85 -K -O >> $PS # original color -I+nt.5
# gmt6 grdimage $topogrd $J $Rmap -Bx -By -Cfes.cpt -I+nt.85 -K -O >> $PS # gray scale

gmt6 pscoast $Rmap $J -Bx -By -Na/.05p -A10 -P -K -Di -O -W.1p >> $PS #-A10+l -Ia
########
# gmt6 psxy $GAB $Rmap $J -W0.3p,navy -Gdarkgrey@70 -O -K >> $PS ###

# gmt6 psxy ~/Research/Lake_eyre_data/geological_ft/CrustalBoundaries.txt $Rmap $J -W1.05p,gray -O -K -V >> $PS ### crustal boundary
# gmt6 psxy ~/Research/Lake_eyre_data/geological_ft/Gawler_Final.gmt $Rmap $J -W0.85p,darkblue,- -O -K >> $PS ### crustal boundary


# awk '{print $1,$2,$3*.1}' GA_18_20_eq.txt | gmt6 psxy $Rmap $J -: -Sc -W.15,black -Gwhite@10 -B -O -P -V -K >> $PS # mine
####
#P_T #darkolivegreen for BW
# darkslategray for color topo
###round_1 fms
# gmt6 psxy ~/Research/Lake_eyre_data/station/old_mix_stuff/stations_all_BB.txt -: -W.05 -St.22 -Gdarkgreen@20 $J $Rmap -O -K >> $PS
# gmt6 psxy ~/Research/Lake_eyre_data/station/old_mix_stuff/stations_all_SP.txt -: -W.05 -St.22 -Gdarkgreen@20 $J $Rmap -O -K >> $PS

awk '{print $1,$2}' 5g_stations_LE.txt | gmt6 psxy -W.01 -Gwhite -St.125  $J $Rmap -O -V -K >> $PS
awk '{print $1,$2}' 6k_stations.txt | gmt6 psxy -W.01 -Gwhite -St.125  $J $Rmap -O -V -K >> $PS


gmt6 psbasemap -B -V $J $Rmap -O -P -K >> $PS

 #-A10+l -Ia
gmt6 gmtset FONT_ANNOT_PRIMARY 8p MAP_FRAME_PEN .8p FONT_LABEL 7.5p
gmt6 psscale -Dx1c/1c+w2.6c/.24c+e+h -O -G-100/880 -F+gwhite+p.1 -Cafrikakarte.cpt -Bx250 -By+l"(m)" >> $PS #-G-100/1000


####LEGEND
# gmt6 pslegend -DJLT+w2.15c/.27c+o-2.15c/-.27c -F+gwhite@10+p.05 -O $J $Rmap << EOF >> $PS
# # #S 0.2c s 0.3c DarkOrange - 0.2i 1Ds-1Dp
# EOF

#pstext $J $Rmap -P -O -K <<End >> $PS
#-40 -7 18 0 31 1 (b)
#End

gmt6 ps2raster -A -Tf -E720 -P -Z -Vq $PS
open map_oz.pdf
# cp map.pdf ~/Dropbox/map_na_aust_ml1.pdf
###
