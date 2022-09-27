#!/bin/csh
# this plots basement depth from drillhole data with Rf estimated basement depth

set PS=coverTh_rf.ps
set Rmap=-R127.5/142/-38.5/-24.5 # SA

set topogrd = /Users/shubham/Research/etopo/etopo_SA.grd
set J=-JM12/11 #sa

set grd=temp.grd

gmt6 gmtset  MAP_FRAME_TYPE fancy+ MAP_FRAME_WIDTH 1p FONT_ANNOT_PRIMARY 8.5p MAP_TICK_LENGTH_PRIMARY 2.5p MAP_FRAME_PEN 0.8p

gmt6 psbasemap -BneWS -Bxa2f1 -Bya2f1 $J $Rmap -X2c -Y2c -K >! $PS

set cover = /Users/shubham/Research/Lake_eyre_data/geological_ft/drillholes_depthtobasement_shp/cover_thickness.txt
set cover = /Users/Shubham/Research/geological_data_straya/GSSA_basement/cover_thickness_2020.txt
set province_nosedi=/Users/Shubham/Research/geological_data_straya/GeoRegions_caroline/province_noSedi.gmt

gmt6 pscoast -Bx -By -A100 -Na/.05p -Dh -W.5,black -Saliceblue -GWhiteSmoke@30 $J $Rmap -O -K -P -V >> $PS

gmt6 psxy $cover $Rmap $J -Ss.2 -Cbatlow_edi_3000.cpt -O -K -V >> $PS ### cover thickness batlow
gmt6 psxy $province_nosedi $Rmap $J -W0.4p,black,- -O -K >> $PS ### crustal boundary
awk '{print $1,$2}' P-delay/stations/*.txt | gmt6 psxy $Rmap $J -W.1 -Sc.04 -Gblack -O -P -K -V >> $PS
gmt6 pslegend -Dx2c/1.75c+w4.2c/.6c+o-1c/-.5c -F+gwhite -O -P -K -V $J $Rmap << EOF >> $PS
S 0.2c s 0.25c pink@10 - 0.2i Borehole basement depth
EOF
awk '{if ($3<0.575  && $3!="NaN") print($1,$2,$3*.37*1000) }' P-delay/P_delay_SA_all.txt | gmt6 psxy  -: -St.38 -W.01,black -Cbatlow_edi_3000.cpt $J $Rmap -O -K -V >> $PS
awk '{if ($3>0.575  && $3!="NaN") print($1,$2,$3*3210-1660) }' P-delay/P_delay_SA_all.txt | gmt6 psxy  -: -St.38 -W.01,black -Cbatlow_edi_3000.cpt $J $Rmap -O -V >> $PS

gmt6 psconvert -A -Tf -P -Z $PS
open coverTh_rf.pdf
