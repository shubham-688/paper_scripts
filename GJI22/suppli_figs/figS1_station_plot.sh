#!/bin/csh

set PS=map.ps

set Rmap=-R127.5/142/-38.5/-24.5 # SA

set J=-JM12/11 #sa

set topogrd = /Users/shubham/Research/etopo/etopo_SA.grd

set topoxyz = ~/Research/etopo/etopo1.xyz
gmt5 gmtset  MAP_FRAME_TYPE fancy+ MAP_FRAME_WIDTH 1p FONT_ANNOT_PRIMARY 8.5p MAP_TICK_LENGTH_PRIMARY 2.5p MAP_FRAME_PEN 0.8p

gmt5 psbasemap -BneWS -Bxa2f1 -Bya2f1 -V $J $Rmap -K >! $PS

set fes = ~/Documents/ScientificColourMaps7/fes/fes.cpt

gmt5 grdimage $topogrd $J $Rmap -Bx -By -Cfes.cpt -I+nt.85 -K -O >> $PS # -I+nt.6 original..increased for extra contarst
gmt5 pscoast $Rmap $J -Bx -By -Na/.05p -A10 -P -K -Di -O -W.1p >> $PS #-A10+l -Ia

gmt5 psxy ~/Research/south_aus/AU_stations_SA.txt -Ss.3 -Gblack $J $Rmap -O -K >> $PS
awk '{print $1-.32,$2,$3}' ~/Research/south_aus/AU_stations_SA.txt | gmt5 pstext $Rmap $J -F+f2.5p,Helvetica-Bold -Gwhite -O -P -K >> $PS
awk '{print $1,$2}' Marla/Marla_stations_SA.txt | gmt5 psxy -W.05 -GDimGray -Si.12  $J $Rmap -O -V -K >> $PS
awk '{print $1,$2}' ASR/ASR_stations_SA.txt | gmt5 psxy -W.1 -GLavender -St.3  $J $Rmap -O -V -K >> $PS
awk '{print $1,$2-.15,$3}' ASR/ASR_stations_SA.txt | gmt5 pstext $Rmap $J -F+f2.8p,Helvetica-Bold -Gwhite -O -P -K >> $PS
awk '{print $1,$2}' Bilby/Bilby_stations_SA.txt | gmt5 psxy -W.1 -GFireBrick -St.3  $J $Rmap -O -V -K >> $PS
awk '{print $1,$2-.15,$3}' Bilby/Bilby_stations_SA.txt | gmt5 pstext $Rmap $J -F+f2.8p,Helvetica-Bold -Gwhite -O -P -K >> $PS
awk '{print $1,$2}' Skippy/Skippy_stations_SA.txt | gmt5 psxy -W.1 -Gcoral -St.3  $J $Rmap -O -V -K >> $PS
awk '{print $1,$2-.15,$3}' Skippy/Skippy_stations_SA.txt | gmt5 pstext $Rmap $J -F+f2.8p,Helvetica-Bold -Gwhite -O -P -K >> $PS
awk '{print $1,$2}' SOC_7K/SOC7K_stations_SA.txt | gmt5 psxy -W.1 -Gindianred -St.3  $J $Rmap -O -V -K >> $PS
awk '{print $1,$2-.15,$3}' SOC_7K/SOC7K_stations_SA.txt | gmt5 pstext $Rmap $J -F+f2.8p,Helvetica-Bold -Gwhite -O -P -K >> $PS
awk '{print $1,$2}' Gawler/Gawler_stations_SA.txt | gmt5 psxy -W.1 -GLightPink -St.3  $J $Rmap -O -V -K >> $PS
awk '{print $1,$2-.15,$3}' Gawler/Gawler_stations_SA.txt | gmt5 pstext $Rmap $J -F+f2.8p,Helvetica-Bold -Gwhite -O -P -K >> $PS
awk '{print $1,$2}' 1F_curna/1F_stations_SA.txt | gmt5 psxy -W.1 -Gtan -St.3  $J $Rmap -O -V -K >> $PS
awk '{print $1,$2-.15,$3}' 1F_curna/1F_stations_SA.txt | gmt5 pstext $Rmap $J -F+f2.8p,Helvetica-Bold -Gwhite -O -P -K >> $PS
gmt5 psxy 5g_stations_LE.txt -W.1 -St.3 -W.1 -GOliveDrab $J $Rmap -O -K >> $PS #
awk '{print $1,$2-.15,$3}' 5g_stations_LE.txt | gmt5 pstext $Rmap $J -F+f2.8p,Helvetica-Bold -Gwhite -O -P -K >> $PS
awk '{print $1,$2}' YJ_capral/YJ_stations_SA.txt | gmt5 psxy -W.1 -GMediumPurple -St.3  $J $Rmap -O -V -K >> $PS
awk '{print $1,$2-.15,$3}' YJ_capral/YJ_stations_SA.txt | gmt5 pstext $Rmap $J -F+f2.8p,Helvetica-Bold -Gwhite -O -P -K >> $PS
awk '{print $1,$2}' Ausis/AUSIS_stations_SA.txt | gmt5 psxy -W.1 -GKhaki -St.3  $J $Rmap -O -V -K >> $PS
awk '{print $1,$2-.15,$3}' Ausis/AUSIS_stations_SA.txt | gmt5 pstext $Rmap $J -F+f2.8p,Helvetica-Bold -Gwhite -O -P -K >> $PS
awk '{print $1,$2}' AQT_1Q/AQT1Q_stations_SA.txt | gmt5 psxy -W.1 -GSandyBrown -St.3  $J $Rmap -O -V -K >> $PS
awk '{print $1-.07,$2-.15,$3}' AQT_1Q/AQT1Q_stations_SA.txt | gmt5 pstext $Rmap $J -F+f2.8p,Helvetica-Bold -Gwhite -O -P -K >> $PS
awk '{print $1,$2}' SA_array/6k_stations.txt | gmt5 psxy -W.1 -GSteelBlue -St.3  $J $Rmap -O -V -K >> $PS
awk '{print $1-.07,$2-.15,$3}' SA_array/6k_stations.txt | gmt5 pstext $Rmap $J -F+f2.5p,Helvetica-Bold -Gwhite -O -P -K >> $PS
awk '{print $1,$2}' tasmal/7I_stations.txt | gmt5 psxy -W.1 -Gdarkcyan -St.3  $J $Rmap -O -V -K >> $PS
awk '{print $1-.07,$2-.15,$3}' tasmal/7I_stations.txt | gmt5 pstext $Rmap $J -F+f2.5p,Helvetica-Bold -Gwhite -O -P -K >> $PS
#####
gmt5 pslegend -Dx1.6c/1.1c+w3.1c/4.85c+o-1c/-.5c -F+gwhite+p.5 -O $J $Rmap << EOF >> $PS
S 0.2c t 0.3c OliveDrab - 0.2i Lake Eyre (5G)
S 0.2c t 0.3c SteelBlue - 0.2i AusArray SA (6K)
S 0.2c t 0.3c Lavender - 0.2i ASR (5J)
S 0.2c t 0.3c coral - 0.2i Skippy (7B)
S 0.2c t 0.3c tan - 0.2i Curnamona (1F)
S 0.2c t 0.3c FireBrick - 0.2i Bilby (6F)
S 0.2c t 0.3c MediumPurple - 0.2i Capral (YJ)
S 0.2c t 0.3c Khaki - 0.2i AUSIS (S1)
S 0.2c t 0.3c indianred - 0.2i SOC (7K)
S 0.2c t 0.3c darkcyan - 0.2i Tasmal (7I)
S 0.2c t 0.3c LightPink - 0.2i Gawler (1G)
S 0.2c t 0.3c SandyBrown - 0.2i AQT (1Q)
S 0.2c i 0.25c DimGray - 0.2i Marla (3G)
S 0.2c s 0.3c Black - 0.2i ANSN (AU)
# S 0.2c c 0.22c firebrick - 0.2i Eq
EOF

gmt5 psconvert -A -Tf -P -Z $PS
open map.pdf
