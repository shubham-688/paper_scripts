#!/bin/csh

set PS=Pdelay.ps
set Rmap_=-R112/157.5/-44/-10
set J_=-JM15
set topogrd = /Users/shubham/Research/etopo/etopo_SA.grd

set Rmap=-R127.5/142/-38.5/-24.5 # SA
set J=-JM12/11 #sa
set grd=temp.grd

set ccpt = ~/Documents/ScientificColourMaps5/cork/cork.cpt
set lake_eyre = ~/Research/geological_data_straya/hydrological_data/Vector_data/lake_eyre.gmt
set province=/Users/Shubham/Research/geological_data_straya/GeoRegions_caroline/shapefiles/province.gmt
set province_nosedi=/Users/Shubham/Research/geological_data_straya/GeoRegions_caroline/province_noSedi.gmt

gmt6 gmtset  MAP_FRAME_TYPE fancy+ MAP_FRAME_WIDTH 1p FONT_ANNOT_PRIMARY 8.5p MAP_TICK_LENGTH_PRIMARY 2.5p MAP_FRAME_PEN 0.8p
gmt6 makecpt -Cabyss.cpt -T-300/0/25 > abyss1.cpt

set oleron_abyss = ~/Research/south_aus/oleron_abyss.cpt

gmt6 psbasemap -BneWS -Bxa2f1 -Bya2f1 $J $Rmap -X2c -Y2c -K >! $PS
####contour
gmt6 grdimage $topogrd $J $Rmap -Bx -By -C$oleron_abyss -I+nt.85 -K -O >> $PS # -I+nt.6 original..increased for extra contarst

gmt6 pscoast -Bx -By -A100 -Na/.05p -Dh -W.5,black $J $Rmap -O -K -P -V >> $PS
####################################
gmt6 psxy $province_nosedi $Rmap $J -W0.4p,DimGray -GLightGrey@50 -O -K >> $PS ### crustal boundary

gmt6 psxy P_delay*.txt -: -W.2 -St.35 -W.1,black -Ccork.cpt -Gred $J $Rmap -O -K >> $PS
gmt6 psxy stations/*.txt -W.2 -St.35 -W.1,black $J $Rmap -O -K >> $PS #
gmt6 psxy $lake_eyre $Rmap $J -W0.02p,dimgray -O -K  >> $PS # lake eyre boundary

awk '{print $1,$2-.17,$3}' stations/*.txt | gmt6 pstext $Rmap $J -F+f2.5p,Helvetica-Bold -Gwhite -O -P -K >> $PS
awk '{print $1,$2-.17,$3}' stations/stations_imp.txt | gmt6 pstext $Rmap $J -F+f2.5p,Helvetica-Bold -Gforestgreen@50 -O -P -K >> $PS

echo 137.65 -28.35 KT-LE | gmt6 pstext $Rmap $J -F+f7p,Courier-bold -O -P -K  >> $PS

gmt6 pslegend -Dx1.5c/1c+w5.2c/1.2c+o-1c/-.5c -F+gwhite+p.7 -O -K $J $Rmap << EOF >> $PS
EOF
gmt6 psscale -Dx.84c/1.05c+w3.c/.25c+e+h -O -Ccork.cpt -Bx.3 -By+l"TPsb (s)" -P >> $PS

gmt6 psconvert -A -Tf -P -Z $PS
open Pdelay.pdf
