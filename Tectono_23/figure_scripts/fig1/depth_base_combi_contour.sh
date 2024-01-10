#!/bin/csh
# this plots basement depth from drillhole data
# awk 'NR==FNR{a[$4]=$0; next} {if ($1 in a){print $1,$2,a[$1]}}' station_try.dat X.txt | awk '{print($3,$4,$2,$1)}' > X_.txt
# set filename='~/Research/AuSREM-CM-moho/AusMoho2012.xyz'

set PS=base_thick_comb.ps
set Rmap=-R130/142/-33/-24 # Lake_eyre+ Flinders
set Rmap=-R130/142/-34/-24 # SA

# set J_=-JM15
set topogrd = /Users/shubham/Research/etopo/etopo_sa_15s.grd
# set topogrd = /Users/shubham/Research/etopo/etopo_SA.grd

set stations_6k = ~/Research/Lake_eyre_data/sa_array_rf/maps/6k_stations.txt
set province_nosedi=/Users/Shubham/Research/geological_data_straya/GeoRegions_caroline/province_noSedi.gmt
set GAB = ~/Research/geological_data_straya/artesian_basin/81672/GAB_Hydrological_Boundary.gmt

#set J=-JM8/22
set J=-JM12/11 #sa

set grd=temp.grd

# gmt6 gmtset  MAP_FRAME_TYPE plain MAP_FRAME_WIDTH 1.5p FONT_ANNOT_PRIMARY 8p MAP_TICK_LENGTH_PRIMARY 1.5p MAP_FRAME_PEN 0.7p FONT_ANNOT_SECONDARY 7p
gmt6 gmtset MAP_FRAME_TYPE plain MAP_FRAME_WIDTH 2p FONT_ANNOT_PRIMARY 8.5p MAP_TICK_LENGTH_PRIMARY 3.5p MAP_FRAME_PEN 1.3p

# gmt6 makecpt -T30/55/2.5 -I -D -Ccopper > $cpt

##########
# cork used for P P_delay
# abyss for bathy

# gmt6 makecpt -C$lajo -T-40/2000/50 -D > cork1.cpt

gmt6 psbasemap -BnEwS -Bxa2f1 -Bya2f1 $J $Rmap  -X2c -Y2c -K >! $PS

set cover = /Users/shubham/Research/Lake_eyre_data/geological_ft/drillholes_depthtobasement_shp/cover_thickness.txt

#### following lines for nearneighbor
# gmt6 nearneighbor base_thick_combi_bmean.txt $Rmap -I0.04 -G$grd -S.3d -N3
# nccopy -k 4 $grd grd.grd
# gmt6 grdimage grd.grd $J $Rmap -Bx -By -Cbatlow_edi.cpt -K -O >> $PS
##############
set RRmap=-R129/142/-38.5/-24 # SA

# gmt6 blockmean base_thick_combi_Phanero.txt $Rmap -C -I.05/.05 > base_thick_combi_bmean.txt

############ for surface
# gmt6 surface base_thick_combi_bmean.txt $RRmap -Gsurf.grd -I.05 -T0.25 -V
# nccopy -k 4 surf.grd surf1.grd

gmt6 grdimage surf1.grd $J $Rmap -Bx -By -Cbatlow_edi_3000.cpt -K -O >> $PS
gmt6 grdcontour surf1.grd $J $Rmap -C1000 -A2000+f5+um -Wa0.15p,white -L200 -Wc0.08p,white -O -P -K >> $PS #-Wa5/25 -Wc2/100 -S20 smooth factor


##############
######### spline

# gmt6 greenspline base_thick_combi_bmean.txt $RRmap -I.1 -D1 -Sq -V -Gspline.grd
# gmt6 grdimage spline.grd $J $Rmap -Bx -By -Cbatlow_edi.cpt -K -O >> $PS
#########################

# gmt6 grdimage $topogrd $J $Rmap -Bx -By -Cabyss1.cpt -I+nt.8 -K -O >> $PS

# pscoast $Rmap $J -Ba5f5/a5f5neWS -N1/.005p -A10000 -P -Di -O -W.01p -K >> $PS

# gmt6 psxy ~/Research/Lake_eyre_data/geological_ft/CrustalBoundaries.txt $Rmap $J -W0.45p,black,. -V -O -K >> $PS ### crustal boundary

gmt6 pscoast $Rmap $J -Bx -By -Na/.05p -A10 -P -SWhiteSmoke -GWhiteSmoke@85 -Di -O -K -W.1p >> $PS #-A10+l -Ia

##
# gmt6 psxy ~/Research/Lake_eyre_data/geological_ft/Gawler_Final.gmt $Rmap $J -W0.5p,darkblue@50,- -V -O -K >> $PS ### crustal boundary


# gmt6 psbasemap -BneWS -Bxa2f1 -Bya2f1 $J $Rmap  -X2c -Y2c -K >! $PS

# gmt6 pscoast -Bx -By -A100 -Na/.05p -Dh  -W.5,black -Saliceblue $J $Rmap -O -K -P -V >> $PS
gmt6 psxy $province_nosedi $Rmap $J -W0.3p,black@20,- -O -K >> $PS ### crustal boundary
# gmt6 psxy meso_protereo_carli.txt $Rmap $J -W0.3p,darkred -O -K >> $PS ###

###stationSa
# gmt6 psxy ~/Research/Lake_eyre_data/station/stations_all_BB.txt -: -W.18,dimgray@10 -Sc.25 $J $Rmap -O -K >> $PS
# gmt6 psxy ~/Research/Lake_eyre_data/station/stations_all_SP.txt -: -W.18,white@10 -Sc.25 $J $Rmap -O -K >> $PS
# gmt6 psxy $stations_6k -W.18,white@10 -Sc.25 $J $Rmap -O -K >> $PS
#
awk '{print $1,$2}' 5g_stations_LE.txt | gmt6 psxy -W.45,black@15 -St.3  $J $Rmap -O -V -K >> $PS
awk '{print $1,$2-.15,$3}' 5g_stations_LE.txt | gmt6 pstext $Rmap $J -F+f2.5p,Helvetica-Bold -Gwhite@20 -O -P -K >> $PS

awk '{print $1,$2}' 6k_stations.txt | gmt6 psxy -W.45,black@15 -St.3  $J $Rmap -O -V -K >> $PS
awk '{print $1,$2-.15,$3}' 6k_stations.txt | gmt6 pstext $Rmap $J -F+f2.5p,Helvetica-Bold -Gwhite@20 -O -P -K >> $PS

##### for plotting point data
# gmt6 psxy cover_thick_noCarrie.txt $Rmap $J -Ss.1 -Cbatlow_edi_3000.cpt -O -K -V >> $PS ### cover thickness batlow
# awk '{print $4,$3,$7}' TableS1_gji.txt | gmt6 psxy  $Rmap $J -Ss.1 -Cbatlow_edi_3000.cpt -O -K -V >> $PS ### cover thickness batlow
######


gmt6 gmtset FONT_ANNOT_PRIMARY 5.5p MAP_TICK_LENGTH_PRIMARY 1.5p MAP_FRAME_PEN .6p

# gmt6 psscale -Dx7.15c/-1.2c+w3.c/.25c+e+h -O -K -Ccork.cpt -Bx.3 -By+l"P delay (s)" -P >> $PS
gmt6 psscale -Dx.5c/1.2c+w1.5c/.15c+h+m -O -Cbatlow_edi_shallow.cpt -Bx50 -By+l"(m)"  -V -K  -P >> $PS
gmt6 psscale -Dx.5c/0.7c+w2.c/.15c+e+h  -O -Cbatlow_edi_3000.cpt -Bx1000 -By+l"(m)" -V  -K -P >> $PS
# gmt6 psscale -Dx1.15c/0.7c+w4.c/.25c+e+h  -O -Ccork1.cpt -Bx500 -By+l"Basement Depth (m)" -V -P >> $PS

# echo 138.8 -32.3 ARC | gmt6 pstext $Rmap $J -F+f6.8p,Helvetica -O -P -K  >> $PS
# echo 131 -26.8 MP | gmt6 pstext $Rmap $J -F+f6.8p,Helvetica -O -P -K  >> $PS
# echo 134.8 -31.6 GC | gmt6 pstext $Rmap $J -F+f6.8p,Helvetica -O -P -K  >> $PS
# echo 136.2 -28.3 DPI | gmt6 pstext $Rmap $J -F+f6p,Courier-bold,black+a-65 -O -P -K  >> $PS

echo 131.4 -31.7 Phanerozoic | gmt6 pstext $Rmap $J -F+f8.7p,Helvetica-bold,black -Gwhite@30 -O -P -K  >> $PS
echo 132.2 -32.05 Sedimentary Thickness | gmt6 pstext $Rmap $J -F+f8.7p,Helvetica-bold,black -Gwhite@30 -O -P -K  >> $PS

gmt6 psxy $Rmap $J -W0.35p,black,- -K -V -O << EOF >> $PS
130.55 -32.95
130.55 -33.2
EOF
gmt6 psxy $Rmap $J -W0.35p,black,- -V -O << EOF >> $PS
130.73 -33.2
132.15 -32.95
EOF

gmt6 ps2raster -A -Tj -E920 -P -Z -Vq $PS

open base_thick_comb.jpg
