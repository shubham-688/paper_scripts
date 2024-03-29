#!/bin/csh

set PS=map.ps
#set Rmap=-Rg
#set Rmap=-R-135/40/-90/90
set Rmap=-R131/142/-32/-25 # Lake_eyre final
set Rmap=-R130/142/-34/-24# SA
#set J=-JM8/22
set J=-JM12/11 #sa
#set J=-Jx0.023i/0.029i

set topogrd = /Users/shubham/Research/etopo/etopo_sa_15s.grd
#gmt6 grdcut ETOPO1_Ice_g_gmt6.grd -Getopo_sa.grd -R126/144/-40/-22
set tect_arch = ~/Research/geological_data_straya/tectonic_provinces/archean_gawler.gmt
set tect_prot = ~/Research/geological_data_straya/tectonic_provinces/proterozoic.gmt
set tect_curna = ~/Research/geological_data_straya/tectonic_provinces/curnamona.gmt
set tect_musgrave = ~/Research/geological_data_straya/tectonic_provinces/musgrave.gmt
set tect_inl = ~/Research/Lake_eyre_data/geological_ft/Tectonic_Zones_Inlier.gmt
set lake_eyre = ~/Research/geological_data_straya/hydrological_data/Vector_data/lake_eyre.gmt
set lake_torrens = ~/Research/geological_data_straya/hydrological_data/Vector_data/lake_torrens.gmt
set GAB = ~/Research/geological_data_straya/artesian_basin/81672/GAB_Hydrological_Boundary.gmt
set warb = ~/Research/south_aus/plots_mixed/sedi_basins/basins/early_paleozoic_warburton.txt

# set topogrd = ~/Research/etopo/etopo1.grd
set topoxyz = ~/Research/etopo/etopo1.xyz

# set fes = ~/Documents/ScientificColourMaps7/fes/fes.cpt

# gmt6 gmtset  MAP_FRAME_TYPE fancy+ MAP_FRAME_WIDTH 1p FONT_ANNOT_PRIMARY 6.5p MAP_TICK_LENGTH_PRIMARY .5p MAP_FRAME_PEN 0.8p
gmt6 gmtset  MAP_FRAME_TYPE fancy+ MAP_FRAME_WIDTH 1p FONT_ANNOT_PRIMARY 8.5p MAP_TICK_LENGTH_PRIMARY 2.5p MAP_FRAME_PEN 0.8p

gmt6 psbasemap -BneWS -Bxa2f1 -Bya2f1 -V $J $Rmap -K >! $PS
#gmt6 pscoast  $Rmap $J -B -Na/.005p -Ia -A1000 -P -Sazure -Glightgrey -Di -O -W.01p -K >> $PS

#pscoast -Rg -J -B15g15 -Dc -A10000 -Glightgrey -P -O -W.01p -K >> map.ps
#gmt6 xyz2grd $Rmap $J -I.017 $topoxyz -Gtopo.grd -V

# set devon = ~/Documents/ScientificColourMaps5/devon/devon.cpt
set fes = ~/Documents/ScientificColourMaps7/fes/fes.cpt

# gmt6 makecpt -C$fes -T-800/800/100 > fes.cpt

#psxy stationSa.txt -Sa.52 -h0 -W0.3+cf -Cabc.cpt $J $Rmap -O -V -K >> $PS
# for -Coleron_abyss.cpt, manually removed +ve values from 'gmt6 makecpt -Coleron.cpt -T-7000/7000/100 > oleron_abyss.cpt'

gmt6 grdimage $topogrd $J $Rmap -Bx -By -Cfes.cpt -I+nt.85 -K -O >> $PS # -I+nt.6 original..increased for extra contarst


gmt6 pscoast $Rmap $J -Bx -By -Na/.25p -A10 -P -K -Di -O -W.1p >> $PS #-A10+l -Ia

# gmt6 psxy ~/Research/Lake_eyre_data/geological_ft/CrustalBoundaries.txt $Rmap $J -W0.65p,black,. -O -K >> $PS ### crustal boundary
gmt6 psxy ~/Research/Lake_eyre_data/geological_ft/Gawler_Final.gmt $Rmap $J -W0.1p,lightsalmon -GMistyRose@30 -O -K >> $PS ### crustal boundary

# gmt6 grdcontour moho.grd $J $Rmap -C2.5 -A2.5+f6p+ukm -S15 -Wa.5 -Wc.25 -T -V -O -P -K >> $PS # -S smooth factor
# gmt6 psxy $tect_arch $Rmap $J -W0.25p,firebrick -GMistyRose@30 -V -O -K  >> $PS # archean
gmt6 psxy $tect_musgrave $Rmap $J -W0.25p,DarkSeaGreen -GHoneydew@30 -O -K  >> $PS # proterozoic musgrave
gmt6 psxy $tect_curna $Rmap $J -W0.25p,slateblue -Gmediumpurple@80 -O -K  >> $PS # proterozoic curna
# gmt6 psxy $warb $Rmap $J -W0.5p,SteelBlue -GSteelBlue@70 -O -K  >> $PS
gmt6 psxy $GAB $Rmap $J -W0.3p,Goldenrod -GLightGoldenrodYellow@60 -O -K >> $PS ###
gmt6 psxy $tect_inl $Rmap $J -W0.25p,peru -Gperu@30 -O -K >> $PS ### crustal boundary

# gmt6 psxy $lake_eyre $Rmap $J -W0.12p,dimgray -O -K  >> $PS # lake eyre boundary
gmt6 psxy $lake_torrens $Rmap $J -W0.12p,dimgray -O -K  >> $PS # lake_torrens boundary

awk '{print $3,$2,$4}' ~/Research/AuSREM-CM-moho/H-Mfiles/EC_GA.xyz | gmt6 psxy -W1.0p,darkred $J $Rmap -O -K >> $PS
awk '{print $3,$2,$4}' ~/Research/AuSREM-CM-moho/H-Mfiles/CF2.xyz | gmt6 psxy -W1.0p,darkred $J $Rmap -O -K >> $PS
awk '{print $3,$2,$4}' ~/Research/AuSREM-CM-moho/H-Mfiles/OD1.xyz | gmt6 psxy -W1.0p,darkred $J $Rmap -O -K >> $PS
awk '{print $3,$2,$4}' ~/Research/AuSREM-CM-moho/H-Mfiles/OD2.xyz | gmt6 psxy -W1.0p,darkred $J $Rmap -O -K >> $PS
awk '{print $3,$2,$4}' ~/Research/AuSREM-CM-moho/H-Mfiles/curna_1.txt | gmt6 psxy -W1.0p,darkred $J $Rmap -O -K >> $PS
awk '{print $3,$2,$4}' ~/Research/AuSREM-CM-moho/H-Mfiles/curna_act.txt | gmt6 psxy -W1.0p,darkred $J $Rmap -O -K >> $PS

gmt6 psxy birds.txt -W0.75p,grey $J $Rmap -P -O -K >> $PS
echo 139.75 -28 Birdsville | gmt6 pstext $Rmap $J -F+f5.5p,Helvetica,dimgray -O -P -K  >> $PS
echo 139.75 -28.2 Track | gmt6 pstext $Rmap $J -F+f5.5p,Helvetica,dimgray -O -P -K  >> $PS

gmt6 psxy $J $Rmap -Sv0.1i+ea -P -O -K -Ggrey -W0.25p << EOF >> $PS
139.35 -28 150 .25i
EOF
# awk '{print $3,$2,$4}' ~/Research/AuSREM-CM-moho/H-Mfiles/*refract*.xyz | gmt6 psxy -W0.25p,maroon $J $Rmap -O -K >> $PS


#echo 141.2 -36 VIC | gmt6 pstext $Rmap $J -F+f8.5p,Helvetica-Bold,MidnightBlue+a90 -O -P -K  >> $PS
echo 141.2 -32.8 NSW | gmt6 pstext $Rmap $J -F+f8.5p,Helvetica-Bold,MidnightBlue+a90 -O -P -K  >> $PS
echo 141.2 -26.5 QLD | gmt6 pstext $Rmap $J -F+f8.5p,Helvetica-Bold,MidnightBlue+a90 -O -P -K  >> $PS
echo 131.5 -25.8 NT | gmt6 pstext $Rmap $J -F+f8.5p,Helvetica-Bold,MidnightBlue -O -P -K  >> $PS
echo 128.8 -28 WA | gmt6 pstext $Rmap $J -F+f8.5p,Helvetica-Bold,MidnightBlue+a90 -O -P -K  >> $PS

echo 133.25 -31.5 Gawler Craton | gmt6 pstext $Rmap $J -F+f8.5p,Helvetica,firebrick -O -P -K  >> $PS
echo 132 -26.5 Musgrave Block | gmt6 pstext $Rmap $J -F+f8.5p,Helvetica,DarkGreen -O -P -K  >> $PS
echo 140.8 -31.5 Curnamona Block | gmt6 pstext $Rmap $J -F+f7p,Helvetica,slateblue -O -P -K  >> $PS
echo 138.6 -32 Flinders Ranges | gmt6 pstext $Rmap $J -F+f6.5p,Helvetica+a-80 -Gwhite@40 -O -P -K  >> $PS
echo 140.5 -28.9 Eromanga Basin | gmt6 pstext $Rmap $J -F+f8.5p,Helvetica,Goldenrod -O -P -K  >> $PS
# echo 138.5 -26.3 Warburton B. | gmt6 pstext $Rmap $J -F+f6.5p,Palatino-Bold,SteelBlue -O -P -K  >> $PS
echo 137 -25.1 Simpson Desert | gmt6 pstext $Rmap $J -F+f6.5p,Palatino-Bold,black -O -P -K  >> $PS

###########
# gmt6 psxy $Rmap $J -W.5p,black -V -O -K << EOF >> $PS
# 133 -26
# 140 -26
# EOF
# echo 132.75 -26 A-A\' | gmt6 pstext $Rmap $J -F+f8.5p,Helvetica,black -O -P -K  >> $PS
# # echo 140 -25.8 A\' | gmt6 pstext $Rmap $J -F+f8.5p,Helvetica -O -P -K  >> $PS
# ####
# gmt6 psxy $Rmap $J -W.5p,black -V -O -K << EOF >> $PS
# 139.8 -26
# 138 -30
# EOF
# echo 140 -26.2 B-B\' | gmt6 pstext $Rmap $J -F+f8.5p,Helvetica,black+a60 -O -P -K  >> $PS
# # echo 138.2 -30 B\' | gmt6 pstext $Rmap $J -F+f8.5p,Helvetica -O -P -K  >> $PS
# ####
# gmt6 psxy $Rmap $J -W.5p,black -V -O -K << EOF >> $PS
# 134.74 -26
# 138 -30
# EOF
# echo 134.6 -25.6 C1-C1\' | gmt6 pstext $Rmap $J -F+f8.5p,Helvetica,black+a-60 -O -P -K  >> $PS
####
gmt6 psxy $Rmap $J -W.45p,navy@15,- -V -O -K << EOF >> $PS
133.8 -25.5
138.5 -31.75
EOF
echo 133.65 -25.3 C-C\' | gmt6 pstext $Rmap $J -F+f7.5p,Helvetica,navy+a-60 -Gwhite@40 -O -P -K  >> $PS
####
# gmt6 psxy $Rmap $J -W.45p,navy@15,- -V -O -K << EOF >> $PS
# 132.5 -30.5
# 140.4 -26.6
# EOF
# echo 132.4 -30.4 D-D\' | gmt6 pstext $Rmap $J -F+f5.5p,Helvetica,navy -O -P -K  >> $PS

# awk '{print $2,$3}' ~/Research/Lake_eyre_data/station/marla.txt | gmt6 psxy -: -Si.15 -GDimGray $J $Rmap -O -K >> $PS ### Marla
awk '{print $1,$2}' 5g_stations_LE.txt | gmt6 psxy -W.1 -Gkhaki -St.3  $J $Rmap -O -V -K >> $PS
awk '{print $1-.07,$2-.15,$3}' 5g_stations_LE.txt | gmt6 pstext $Rmap $J -F+f2.5p,Helvetica-Bold -Gwhite -O -P -K >> $PS

awk '{print $1,$2}' 6k_stations.txt | gmt6 psxy -W.1 -GOliveDrab -St.3  $J $Rmap -O -V -K >> $PS
awk '{print $1-.07,$2-.15,$3}' 6k_stations.txt | gmt6 pstext $Rmap $J -F+f2.5p,Helvetica-Bold -Gwhite -O -P -K >> $PS


# awk '{print $1,$2}' AU_stations_SA.txt | gmt6 psxy -W.1 -GSteelBlue -St.3  $J $Rmap -O -V -K >> $PS
# awk '{print $1-.07,$2-.15,$3}' AU_stations_SA.txt | gmt6 pstext $Rmap $J -F+f2.5p,Helvetica-Bold -Gwhite -O -P -K >> $PS

# gmt6 psscale -Dx11.2c/1.4c+w3c/.3c -O -K -G-700/700 -Coleron1.cpt -Bx200 -By+l"m" >> $PS #-G-100/1000

echo 137.3 -28.45 KT-LE | gmt6 pstext $Rmap $J -F+f4p,Courier -O -P -K  >> $PS

echo 131 -30.35 EC-GW | gmt6 pstext $Rmap $J -F+f6p,Helvetica-Bold,darkred+a10 -O -P -K  >> $PS
echo 133.85 -28 GOMA | gmt6 pstext $Rmap $J -F+f6p,Helvetica-Bold,darkred+a-75 -O -P -K  >> $PS
echo 139 -24.7 14-CF2 | gmt6 pstext $Rmap $J -F+f6p,Helvetica-Bold,darkred -O -P -K  >> $PS
echo 137.4 -29.45 OD1 | gmt6 pstext $Rmap $J -F+f6p,Helvetica-Bold,darkred -O -P -K  >> $PS
echo 137.2 -30.6 OD2 | gmt6 pstext $Rmap $J -F+f5p,Helvetica-Bold,darkred -O -P -K  >> $PS
echo 140.5 -30.6 09-C1 | gmt6 pstext $Rmap $J -F+f5.5p,Helvetica-Bold,darkred -O -P -K  >> $PS
echo 140 -31.8 03-CU | gmt6 pstext $Rmap $J -F+f5.5p,Helvetica-Bold,darkred -O -P -K  >> $PS

### towns
gmt6 psxy towns.txt -Sc.1 -: -Gblack $J $Rmap -O -P -K >> $PS

gmt6 pstext -: $Rmap $J -F+f3.5p,Helvetica-Bold -O -K -P <<eof >> $PS
-32.35 137.765833 Port Augusta
-28.85 134.755556 Coober Pedy
-30.43333 138.4 Leigh Creek
eof
awk '{print $1,$2+.41,$3}' towns_s.txt | gmt6 pstext -: $Rmap $J -F+f3.5p,Helvetica-Bold -O -P -K >> $PS
##########

gmt6 psimage fig_inset.png -Dx12.2c/11c+w2.75c/2.2c+jTR+w.1i -V -O -K -P >> $PS

gmt6 gmtset FONT_ANNOT_PRIMARY 7.5p MAP_FRAME_PEN .8p FONT_LABEL 7.5p

gmt6 pslegend -Dx1.6c/1.1c+w3.1c/.85c+o-1c/-.5c -F+gwhite+p.5 -O $J $Rmap << EOF >> $PS
S 0.2c t 0.3c Khaki - 0.2i Lake Eyre (5G)
S 0.2c t 0.3c OliveDrab - 0.2i AusArray SA (6K)
#S 0.2c s 0.3c Black - 0.2i ANSN (AU)
# S 0.2c c 0.22c firebrick - 0.2i Eq
EOF

####LEGEND
# gmt6 pslegend -DJLT+w2.7c/.4c+o-2.7c/-.4c -F+gwhite+p.1 -O $J $Rmap << EOF >> $PS
#S 0.2c s 0.3c gold - 0.2i 1Ds-1Dp
# EOF

# echo Oct18-July19 | gmt6 pstext $J $Rmap -F+cTL+f12p,darkred -O -P >> $PS

#pstext $J $Rmap -P -O -K <<End >> $PS
#-40 -7 18 0 31 1 (b)
#End

gmt6 ps2raster -A -Tj -E920 -P -Z -Vq $PS
open map.jpg
