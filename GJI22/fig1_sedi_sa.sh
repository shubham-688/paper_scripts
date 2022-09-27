#!/bin/csh

# awk 'NR==FNR{a[$4]=$0; next} {if ($1 in a){print $1,$2,a[$1]}}' station_try.dat X.txt | awk '{print($3,$4,$2,$1)}' > X_.txt

set PS=sedi_basin.ps
set topogrd = etopo_sa.grd
set Rmap=-R130/142/-33/-24 # Lake_eyre+ Flinders
set Rmap=-R127.5/142/-38.5/-24.5 # SA
set J=-JM12/11 #sa
set grd=temp.grd
set lake_eyre = ~/Research/geological_data_straya/hydrological_data/Vector_data/lake_eyre.gmt
set lake_torrens = ~/Research/geological_data_straya/hydrological_data/Vector_data/lake_torrens.gmt
set tect_inl = ~/Research/Lake_eyre_data/geological_ft/Tectonic_Zones_Inlier.gmt

gmt6 gmtset  MAP_FRAME_TYPE fancy+ MAP_FRAME_WIDTH 1p FONT_ANNOT_PRIMARY 8.5p MAP_TICK_LENGTH_PRIMARY 2.5p MAP_FRAME_PEN 0.8p

gmt6 psbasemap -BneWS -Bxa2f1 -Bya2f1 $J $Rmap -X2c -Y2c -K >! $PS

# gmt6 xyz2grd temp_moho.txt $Rmap -I0.5 -V -G$grd
# gmt6 nearneighbor temp_moho.txt $Rmap -I0.5 -G$grd -S.5d -N4
# nccopy -k 4 $grd moho_05.grd

# gmt6 grdimage moho_05.grd $J $Rmap -Bx -By -C$cpt -K -O -V >> $PS

gmt6 grdimage $topogrd $J $Rmap -Bx -By -Coleron_abyss.cpt -I+nt.85 -K -O >> $PS # -I+nt.6 original..increased for extra contarst
gmt6 pscoast $Rmap $J -Bx -By -Na/.05p -A10 -P -K -Di -O -W.1p >> $PS #-A10+l -Ia
gmt6 psxy Gawler_Final.gmt $Rmap $J -W0.75p,darkblue,- -O -K >> $PS ### crustal boundary

gmt6 psxy basins/meso_protereo_carli.txt $Rmap $J -W0.5p,maroon -Gmaroon@40 -O -K  >> $PS
# or seinna
gmt6 psxy basins/neo_proterozoic_amadeus.txt $Rmap $J -W0.5p,maroon -Gmaroon@40 -O -K  >> $PS
gmt6 psxy basins/early_paleozoic_warburton.txt $Rmap $J -W0.5p,SteelBlue -GSteelBlue@50 -O -K  >> $PS
gmt6 psxy basins/late_paleozoic_pedrika.txt $Rmap $J -W0.5p,SteelBlue -GSteelBlue@50 -O -K  >> $PS
# gmt6 psxy basins/late_paleozoic_pedrika_neg.txt $Rmap $J -W0.5p,SteelBlue -GSteelBlue@100 -O -K  >> $PS
gmt6 psxy basins/paleo_meso_cooper.txt $Rmap $J -W0.5p,DarkGreen -GDarkGreen@50 -O -K  >> $PS
gmt6 psxy basins/mesozoic_euromanga.txt $Rmap $J -W0.5p,SandyBrown -GSandyBrown@50 -O -K  >> $PS
gmt6 psxy basins/meso_ceno_otway.txt $Rmap $J -W0.5p,SandyBrown -GSandyBrown@50 -O -K  >> $PS
gmt6 psxy basins/cenozoic_Lake_eyre.txt $Rmap $J -W0.5p,LightGoldenrodYellow -GLightGoldenrodYellow@60 -O -K  >> $PS
gmt6 psxy $tect_inl $Rmap $J -W0.15p,white -Gwhite -O -K >> $PS ### crustal boundary
gmt6 pscoast -Bx -By -A100 -Na/.05p -Dh  -W.5,black -Saliceblue $J $Rmap -O -K -P >> $PS
gmt6 psxy $Rmap $J mine_sa.txt -Sc.23 -Gdimgray -B -O -P -K >> $PS
# gmt6 psxy $Rmap $J mine_sa.txt -Sx.28 -W2.5,black -B -O -P -K >> $PS
gmt6 psxy $Rmap $J mine_sa.txt -Sx.18 -W1,white -B -O -P -K >> $PS

####################################################################################
echo 141.2 -36 VIC | gmt6 pstext $Rmap $J -F+f8.5p,Helvetica-Bold,MidnightBlue+a90 -O -P -K  >> $PS
echo 141.2 -30.8 NSW | gmt6 pstext $Rmap $J -F+f8.5p,Helvetica-Bold,MidnightBlue+a90 -O -P -K  >> $PS
echo 141.2 -26.5 QLD | gmt6 pstext $Rmap $J -F+f8.5p,Helvetica-Bold,MidnightBlue+a90 -O -P -K  >> $PS
echo 133.5 -25.8 NT | gmt6 pstext $Rmap $J -F+f8.5p,Helvetica-Bold,MidnightBlue -O -P -K  >> $PS
echo 128.8 -28 WA | gmt6 pstext $Rmap $J -F+f8.5p,Helvetica-Bold,MidnightBlue+a90 -O -P -K  >> $PS

echo 131.6 -34. SEDIMENTARY BASINS | gmt6 pstext $Rmap $J -F+f10p,Helvetica -Gwhite -O -P -K  >> $PS
echo 131.6 -34.5 OF | gmt6 pstext $Rmap $J -F+f10p,Helvetica -Gwhite -O -P -K  >> $PS
echo 131.6 -35. SOUTH AUSTRALIA | gmt6 pstext $Rmap $J -F+f10p,Helvetica -Gwhite -O -P -K  >> $PS
#############################
echo 130 -28.5 Officer B. | gmt6 pstext $Rmap $J -F+f6.5p,Palatino-Bold,Maroon -Gwhite@60 -O -P -K  >> $PS
echo 131.5 -25 Amadeus B. | gmt6 pstext $Rmap $J -F+f6.5p,Palatino-Bold,Maroon -Gwhite@60 -O -P -K  >> $PS
echo 137 -31.8 Cr. B. | gmt6 pstext $Rmap $J -F+f7p,Palatino-Bold,Maroon+a-60 -Gwhite@60 -O -P -K  >> $PS
echo 141.1 -28.2 Cooper B. | gmt6 pstext $Rmap $J -F+f6.5p,Palatino-Bold,DarkGreen -Gwhite@60 -O -P -K  >> $PS
echo 136 -25 Pedrika B. | gmt6 pstext $Rmap $J -F+f6.5p,Palatino-Bold,SteelBlue -Gwhite@60 -O -P -K  >> $PS
echo 138.5 -26.3 Warburton B. | gmt6 pstext $Rmap $J -F+f6.5p,Palatino-Bold,SteelBlue -Gwhite@60 -O -P -K  >> $PS
echo 140 -31 Arrowie B. | gmt6 pstext $Rmap $J -F+f6.5p,Palatino-Bold,SteelBlue -Gwhite@50 -O -P -K  >> $PS
echo 133.8 -29.3 Arckaringa B. | gmt6 pstext $Rmap $J -F+f6.5p,Palatino-Bold,SteelBlue -Gwhite@60 -O -P -K  >> $PS
echo 128.7 -31 Bight B. | gmt6 pstext $Rmap $J -F+f6.5p,Palatino-Bold,Chocolate -Gwhite@60 -O -P -K  >> $PS
echo 141.2 -33.6 Berri B. | gmt6 pstext $Rmap $J -F+f6.5p,Palatino-Bold,Chocolate+a25 -Gwhite@60 -O -P -K  >> $PS
echo 140.8 -37.2 Otway B. | gmt6 pstext $Rmap $J -F+f6.5p,Palatino-Bold,Chocolate -Gwhite@60 -O -P -K  >> $PS
echo 140 -29.7 Eromanga B. | gmt6 pstext $Rmap $J -F+f7p,Palatino-Bold,Chocolate -Gwhite@60 -O -P -K  >> $PS
echo 140 -35.5 Murray B. | gmt6 pstext $Rmap $J -F+f6.5p,Palatino-Bold,Goldenrod -Gwhite@60 -O -P -K  >> $PS
echo 133 -31 Eucla B. | gmt6 pstext $Rmap $J -F+f6.5p,Palatino-Bold,Goldenrod+a-30 -Gwhite@60 -O -P -K  >> $PS
echo 137 -27 Lake Eyre B. | gmt6 pstext $Rmap $J -F+f6.5p,Palatino-Bold,Goldenrod+a-30 -Gwhite@60 -O -P -K  >> $PS
echo 137.3 -28.35 KT- | gmt6 pstext $Rmap $J -F+f7p,Courier-bold -O -P -K  >> $PS
echo 137.25 -28.6 LE | gmt6 pstext $Rmap $J -F+f7p,Courier-bold -O -P -K  >> $PS
echo 137.9 -31 LT | gmt6 pstext $Rmap $J -F+f7p,Courier-bold -O -P -K  >> $PS
echo 136.1 -28.2 DPI | gmt6 pstext $Rmap $J -F+f7p,Courier-bold,black+a-65 -O -P -K  >> $PS
echo 138.6 -32.7 Flinders Ranges | gmt6 pstext $Rmap $J -F+f7.5p,Helvetica+a-80 -Gwhite@40 -O -P -K  >> $PS
echo 131 -26.5 Musgrave Province | gmt6 pstext $Rmap $J -F+f8p,Helvetica -Gwhite@40 -O -P -K  >> $PS
echo 134.8 -31.6 Gawler Craton | gmt6 pstext $Rmap $J -F+f10p,Helvetica -Gwhite@40 -O -P -K  >> $PS
# echo 133.5 -30.3 Gawler Craton | gmt6 pstext $Rmap $J -F+f10p,Helvetica -Gwhite@40 -O -P -K  >> $PS

gmt6 pstext $Rmap $J -F+f4p,Helvetica -Gwhite@40 -O -P -K << eof >> $PS
136.873056 -30.65 Olympic Dam
135.583333 -29.95 Prominent Hill
137.492778 -31.45 Carrapateena
139.506639 -30.4 Four Mile
137.151510 -33.25 Middleback Ranges
#139.6 -28.10838 Moomba
140.5 -32.33 White Dam
132 -31.17 Jacinth Ambrosia
139.65 -35.094716 Kanmantoo
eof

##########################
gmt6 psxy $lake_eyre $Rmap $J -W0.02p,dimgray -O -K  >> $PS # lake eyre boundary
gmt6 psxy $lake_torrens $Rmap $J -W0.02p,dimgray -O -K  >> $PS # lake_torrens boundary
gmt6 psxy towns.txt -Sc.1 -: -Gblack $J $Rmap -O -P -K >> $PS

gmt6 pstext -: $Rmap $J -F+f3.5p,Helvetica-Bold -O -K -P <<eof >> $PS
-37.68 140.782778 Mount Gambier
-32.35 137.765833 Port Augusta
-28.85 134.755556 Coober Pedy
-30.43333 138.4 Leigh Creek
eof
awk '{print $1+.15,$2,$3}' towns_s.txt | gmt6 pstext -: $Rmap $J -F+f3.5p,Helvetica-Bold -O -P -K >> $PS
gmt6 pslegend -Dx1.6c/1.c+w3.85c/1.25c+o-1c/-.5c -F+gwhite -K -O $J $Rmap << EOF >> $PS
#H 5p,Helvetica Age of sediments
S 0.2c s 0.3c LightGoldenrodYellow .1p 0.2i Cenozoic
S 0.2c s 0.3c SandyBrown@10 .1p 0.2i Mesozoic
S 0.2c s 0.3c DarkGreen@10 .1p 0.2i Mesozoic-Paleozoic
EOF

gmt6 pslegend -Dx5.45c/1.c+w2.85c/1c+o-1c/-.5c -F+gwhite -O $J $Rmap << EOF >> $PS
S 0.2c s 0.3c SteelBlue@10 .1p 0.2i Paleozoic
S 0.2c s 0.3c maroon@20 .1p 0.2i Proterozoic
EOF

gmt6 psconvert -A -Tf -P -Z -Vq $PS

open sedi_basin.pdf
