#!/bin/csh
# this plots basement depth from rf on top of SA sedi basins
set PS=coverTh_rf_sedi.ps
set Rmap=-R127.5/142/-38.5/-24.5 # SA

set topogrd = /Users/shubham/Research/etopo/etopo_SA.grd
set J=-JM12/11 #sa

set grd=temp.grd

gmt6 gmtset  MAP_FRAME_TYPE fancy+ MAP_FRAME_WIDTH 1p FONT_ANNOT_PRIMARY 8.5p MAP_TICK_LENGTH_PRIMARY 2.5p MAP_FRAME_PEN 0.8p

gmt6 psbasemap -BnEwS -Bxa2f1 -Bya2f1 $J $Rmap -X2c -Y2c -K >! $PS

set cover = /Users/shubham/Research/Lake_eyre_data/geological_ft/drillholes_depthtobasement_shp/cover_thickness.txt
gmt6 pscoast -Bx -By -A100 -Na/.05p -Dh -W.5,black -Saliceblue $J $Rmap -O -K -P >> $PS

#####
gmt6 psxy sedi_basins/basins/meso_protereo_carli.txt $Rmap $J -W0.5p,maroon -Gmaroon@40 -O -K  >> $PS
# or seinna
gmt6 psxy sedi_basins/basins/neo_proterozoic_amadeus.txt $Rmap $J -W0.5p,maroon -Gmaroon@40 -O -K  >> $PS
gmt6 psxy sedi_basins/basins/meso_ceno_otway.txt $Rmap $J -W0.5p,SandyBrown -GSandyBrown@50 -O -K  >> $PS
gmt6 psxy sedi_basins/basins/cenozoic_Lake_eyre.txt $Rmap $J -W0.5p,LightGoldenrodYellow -GLightGoldenrodYellow@60 -O -K  >> $PS
###
gmt6 pscoast -Bx -By -A100 -Na/.05p -Dh -W.5,black -Saliceblue $J $Rmap -O -K -P >> $PS
awk '{print $1,$2}' P-delay/stations/*.txt | gmt6 psxy $Rmap $J -W.1 -Sc.04 -Gblack -O -P -K -V >> $PS
awk '{print $1,$2-.17,$3}' P-delay/stations/*.txt | gmt6 pstext $Rmap $J -F+f1.5p,Helvetica-Bold -Gwhite -O -P -K >> $PS
awk '{if ($3<0.575  && $3!="NaN") print($1,$2,$3*.37*1000) }' P-delay/P_delay_SA_all.txt | gmt6 psxy  -: -St.28 -W.01,black -Cbatlow_edi_3000.cpt $J $Rmap -O -K -V >> $PS
awk '{if ($3>0.575  && $3!="NaN") print($1,$2,$3*3210-1660) }' P-delay/P_delay_SA_all.txt | gmt6 psxy  -: -St.28 -W.01,black -Cbatlow_edi_3000.cpt $J $Rmap -O -K -V >> $PS


gmt6 psxy $Rmap $J -W0.45p,black,- -K -V -O << EOF >> $PS
129 -37.45
129 -36.82
EOF
gmt6 psxy $Rmap $J -W0.45p,black,- -K -V -O << EOF >> $PS
129.38 -37.45
133 -36.78
EOF
#
gmt6 pstext $Rmap $J -F+f10p -V -O -P -K << eof >> $PS
133 -36.4 Age of sediments
eof
# H 8p,Helvetica
gmt6 pslegend -Dx1.6c/1.c+w3.85c/1.25c+o-1c/-.5c -F+gwhite -K -O $J $Rmap << EOF >> $PS
S 0.2c s 0.3c LightGoldenrodYellow .1p 0.2i Cenozoic
S 0.2c s 0.3c SandyBrown@10 .1p 0.2i Mesozoic
S 0.2c s 0.3c DarkGreen@10 .1p 0.2i Mesozoic-Paleozoic
#L 9p,Times-Roman R Smith et al
EOF
gmt6 pslegend -Dx5.45c/1.25c+w2.85c/1c+o-1c/-.5c -F+gwhite -O $J $Rmap << EOF >> $PS
S 0.2c s 0.3c SteelBlue@10 .1p 0.2i Paleozoic
S 0.2c s 0.3c maroon@20 .1p 0.2i Proterozoic
EOF


gmt6 psconvert -A -Tf -P -Z -Vq $PS

open coverTh_rf_sedi.pdf
