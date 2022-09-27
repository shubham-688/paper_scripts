#!/bin/csh
# this plots basement depth from drillhole data

set PS=coverTh_BH.ps
set Rmap=-R128.75/141.25/-38/-25.75 # SA
set J_=-JM15
set topogrd = /Users/shubham/Research/etopo/etopo_SA.grd

set J=-JM12/11 #sa
set grd=temp.grd
gmt6 gmtset  MAP_FRAME_TYPE fancy+ MAP_FRAME_WIDTH 1p FONT_ANNOT_PRIMARY 8.5p MAP_TICK_LENGTH_PRIMARY 2.5p MAP_FRAME_PEN 0.8p
# gmt6 makecpt -T30/55/2.5 -I -D -Ccopper > $cpt
##########
# abyss for bathy

# gmt6 makecpt -C$lajo -T-40/2000/50 -D > cork1.cpt
gmt6 psbasemap -BneWS -Bxa2f1 -Bya2f1 $J $Rmap -X2c -Y2c -K >! $PS

set cover = /Users/Shubham/Research/geological_data_straya/GSSA_basement/cover_thickness_2020.txt
set province=/Users/Shubham/Research/geological_data_straya/GeoRegions_caroline/shapefiles/province.gmt
set province_nosedi=/Users/Shubham/Research/geological_data_straya/GeoRegions_caroline/province_noSedi.gmt


gmt6 pscoast -Bx -By -A100 -Na/.05p -Dh -W.35,black -Saliceblue -GWhiteSmoke@30 $J $Rmap -O -K -P -V >> $PS
gmt6 psxy cover_thickness_2020.txt $Rmap $J -Ss.2 -Cbatlow_edi_3000.cpt -O -K -V >> $PS ### cover thickness batlow

gmt6 psxy $province_nosedi $Rmap $J -W0.4p,black,- -O -K >> $PS ### crustal boundary
gmt6 pscoast -Bx -By -A100 -Na/.05p -Dh -W.35,black $J $Rmap -O -K -P -V >> $PS

echo 138.8 -32.3 ARC | gmt6 pstext $Rmap $J -F+f9p,Helvetica -O -P -K  >> $PS
echo 131 -26.8 MP | gmt6 pstext $Rmap $J -F+f9p,Helvetica -O -P -K  >> $PS
echo 134.8 -31.6 GC | gmt6 pstext $Rmap $J -F+f9p,Helvetica -O -P -K  >> $PS
echo 136.2 -28.3 DPI | gmt6 pstext $Rmap $J -F+f8p,Courier-bold,black+a-65 -O -P -K  >> $PS

gmt6 psxy $Rmap $J -W0.45p,black,- -K -V -O << EOF >> $PS
130.05 -37.1
130.05 -36.52
EOF
gmt6 psxy $Rmap $J -W0.45p,black,- -K -V -O << EOF >> $PS
130.43 -37.1
133.45 -36.54
EOF


gmt6 psscale -Dx1.15c/1.7c+w3.c/.25c+h+m -O -Cbatlow_edi_shallow.cpt -Bx50 -By+l"(m)"  -V -K  -P >> $PS
gmt6 psscale -Dx1.15c/0.7c+w4.c/.25c+e+h  -O -Cbatlow_edi_3000.cpt -Bx1000 -By+l"Basement Depth (m)" -V  -P >> $PS

gmt6 psconvert -A -Tf -P -Z $PS
#open coverTh_BH.pdf
