#!/bin/csh

set PS=fig_grav.ps
#set Rmap=-Rg
#set Rmap=-R-135/40/-90/90
set Rmap=-R131/142/-32/-25 # Lake_eyre final
set Rmap=-R130/142/-34/-24 # SA
#set J=-JM8/22
set J=-JM12/11 #sa
#set J=-Jx0.023i/0.029i

#ogr2ogr -f "GMT" Hudson_bounds.gmt  Hudson_bounds.shp -t_srs WGS84
set topogrd = /Users/shubham/Research/etopo/etopo_sa_10s.grd
#gmt6 grdcut ETOPO1_Ice_g_gmt6.grd -Getopo_sa.grd -R126/144/-40/-22
set tpsfhm = ~/Documents/cpt-city/ma/gray/grayscale02

# set topoxyz = ~/Research/etopo/etopo1.xyz

# set fes = ~/Documents/ScientificColourMaps7/fes/fes.cpt
set eq_all_LE = ~/Research/earthquakes/ga_database/LE_allmag_alltime_FR.txt
set eq_all_SA = ~/Research/earthquakes/ga_database/SA_mag2_4.txt
set stations_6k = ~/Research/Lake_eyre_data/sa_array_rf/maps/6k_stations.txt

set mag_tiff_hi = ~/Research/Lake_eyre_data/gravity_plots/ga_aus_grav_19.tif

set GRAV=~/Research/GravMag/onshore_geodetic_Complete_Bouguer_2016.nc
set GRAV19=~/Research/GravMag/Gravmap2019-grid-grv_cscba.nc

set GAB = ~/Research/geological_data_straya/artesian_basin/81672/GAB_Hydrological_Boundary.gmt
set basement = ~/Research/Lake_eyre_data/geological_ft/4f_Basement_Inlier/Basement_Inlier.gmt
set faults = ~/Research/geological_data_straya/shapefiles_surafce_geology_australia/faults.gmt
set tect_arch = ~/Research/geological_data_straya/tectonic_provinces/archean_gawler.gmt
set tect_prot = ~/Research/geological_data_straya/tectonic_provinces/proterozoic.gmt
set tect_curna = ~/Research/geological_data_straya/tectonic_provinces/curnamona.gmt
set tect_musgrave = ~/Research/geological_data_straya/tectonic_provinces/musgrave.gmt
set sedi_LE = ~/Research/geological_data_straya/sedi_basin/cenozoic_Lake_eyre.txt

set tect_THZ = ~/Research/Lake_eyre_data/geological_ft/Tectonic_Zones_THZ.gmt
set tect_inl = ~/Research/Lake_eyre_data/geological_ft/Tectonic_Zones_Inlier.gmt

set lake_eyre = ~/Research/geological_data_straya/hydrological_data/Vector_data/lake_eyre.gmt
set lake_torrens = ~/Research/geological_data_straya/hydrological_data/Vector_data/lake_torrens.gmt

set sedi_officer = ~/Research/geological_data_straya/sedi_basin/neo_proterozoic_amadeus.txt
###############################################
###############################################

gmt6 gmtset  MAP_FRAME_TYPE plain MAP_FRAME_WIDTH 10p FONT_ANNOT_PRIMARY 8.5p MAP_TICK_LENGTH_PRIMARY 3.5p MAP_FRAME_PEN 1.3p

gmt6 psbasemap -BNeWs -Bxa2f1 -Bya2f1 $J $Rmap -K >! $PS
#gmt6 grdimage $topogrd $J $Rmap -CETOPO1.cpt -K -O -P  >> $PS
#gmt6 pscoast  $Rmap $J -B -Na/.005p -Ia -A1000 -P -Sazure -Glightgrey -Di -O -W.01p -K >> $PS

#pscoast -Rg -J -B15g15 -Dc -A10000 -Glightgrey -P -O -W.01p -K >> map.ps

set grayC = ~/Documents/ScientificColourMaps7/grayC/grayC.cpt

set tokyo = ~/Documents/ScientificColourMaps7/bilbao/bilbao.cpt

# gmt6 makecpt -Fgray -Cetopo1 -V > etopocolor.cpt # etopo1, -A
# gmt6 makecpt -C$bilbao -D -T4/6.5/.2 > bilbao.cpt
# gmt6 makecpt -Croma.cpt -I -T4/6.5/.25 > devon.cpt


# gmt6 makecpt -C$tpsfhm -T-10000/0/1000  > tpsfhm.cpt
# gmt6 makecpt -C$grayC -T-11000/0/500 -I > gray.cpt
##############
gmt6 makecpt -C$tokyo -T-650/250/10 -D  > bil.cpt
gmt6 makecpt -C$grayC -T-650/250/10 -I -D  > gray.cpt
gmt6 makecpt -Croma.cpt -T30/55/2.5 -I > devon_m.cpt

#############

#psxy stationSa.txt -Sa.52 -h0 -W0.3+cf -Cabc.cpt $J $Rmap -O -V -K >> $PS
# for -Coleron_abyss.cpt, manually removed +ve values from 'gmt6 makecpt -Coleron.cpt -T-7000/7000/100 > oleron_abyss.cpt'

# gmt6 grdimage $topogrd $J $Rmap -Bx -By -Ctpsfhm.cpt -I+nt.85 -K -O >> $PS # -I+nt.6 original..increased for extra contarst

# gmt6 grdimage $mag_tiff_hi $J $Rmap -M -Gb -Bx -By -V -K -O >> $PS # here .65 reduce the intensity by 65%.
## -M forces Monochrome
# gmt6 grdgradient $GRAV19 -Ggradient.grd=nb/a -A0/270 -Ne0.5 -V
gmt6 grdimage $GRAV19 $J $Rmap -Cbil.cpt -I+nt.8 -Bx -By -V -K -O >> $PS # -M is for monochrome

# gmt6 psxy ~/Research/Lake_eyre_data/geological_ft/CrustalBoundaries.txt $Rmap $J -W0.65p,black,- -O -K >> $PS ### crustal boundary


# gmt6 psxy $lake_eyre $Rmap $J -W0.09p,dimgray -O -K  >> $PS # lake eyre boundary
# gmt6 psxy $lake_torrens $Rmap $J -W0.09p,dimgray -O -K  >> $PS # lake_torrens boundary

# echo 137.3 -28.4 Lake | gmt6 pstext $Rmap $J -F+f4p,Courier -O -P -K  >> $PS
# echo 137.3 -28.55 Eyre | gmt6 pstext $Rmap $J -F+f4p,Courier -O -P -K  >> $PS

# gmt6 psxy ~/Research/Lake_eyre_data/geological_ft/CrustalBoundaries.txt $Rmap $J -W0.6p,white,- -O -K >> $PS ### crustal boundary

# gmt6 psxy $basement $Rmap $J -W0.15p,darkred -O -K >> $PS ### crustal boundary
# gmt6 psxy $tect_inl $Rmap $J -W0.15p,rosybrown -Grosybrown@30 -O -K >> $PS ### Inlier


# gmt6 psxy ~/Research/south_aus/AU_stations_SA.txt -Ss.2 -Gblack $J $Rmap -O -K >> $PS
gmt6 psxy ~/Research/Lake_eyre_data/station/old_mix_stuff/stations_all_BB.txt -: -W.08 -Sc.25 -Gwhite $J $Rmap -O -K >> $PS
gmt6 psxy ~/Research/Lake_eyre_data/station/old_mix_stuff/stations_all_SP.txt -: -W.08 -Sc.25 -Gwhite $J $Rmap -O -K >> $PS
gmt6 psxy $stations_6k -W.08 -Sc.25 -Gwhite $J $Rmap -O -K >> $PS

# awk '{print $1,$2,($3-$4)*8.9}' stations_all_BB_Pms.txt > temp.txt

gmt6 psxy moho_oneVal.txt -: -W.08 -Sc.32 -Cdevon_m.cpt $J $Rmap -O -K >> $PS # stations with one Moho
awk '{print $1,$2,$4}' stations_all_moho1_2.txt | gmt6 psxy -: -W.08 -Sc.32 -Cdevon_m.cpt $J $Rmap -O -K >> $PS # plots higher pms value
gmt6 psxy stations_all_moho1_2.txt -: -W.08 -Sc.19 -Cdevon_m.cpt $J $Rmap -O -K >> $PS # # plots lower pms value
# awk '{print $2,$1-.22,$3}' ~/Research/Lake_eyre_data/station/old_mix_stuff/stations_all_*.txt | gmt6 pstext $Rmap $J -F+f2.5p,Helvetica-Bold -Gwhite -O -P -K >> $PS

#
gmt6 psxy 6k_stations_Ps-dt.txt -: -W.08 -Sc.32 -Cdevon_m.cpt $J $Rmap -O -K >> $PS # stations with one Moho
awk '{print $1,$2,$4*8.9}' 6k_stations_Ps.txt | gmt6 psxy -W.08 -Cdevon_m.cpt -Sc.32  $J $Rmap -O -V -K >> $PS
# awk '{print $1-.07,$2-.23,$3}' $stations_6k | gmt6 pstext $Rmap $J -F+f2.5p,Helvetica-Bold -Gwhite -O -P -K >> $PS

#########
# Passive
# awk '{print $3,$2,$4}' ~/Research/AuSREM-CM-moho/H-Mfiles/H-HKbk.xyz | gmt6 psxy -W.01 -Ss.18 -Clajolla.cpt $J $Rmap -O -K >> $PS # Not in SA

awk '{print $3,$2,$4}' ~/Research/AuSREM-CM-moho/H-Mfiles/*AC*.xyz | gmt6 psxy -W.001 -Ss.18 -Cdevon_m.cpt $J $Rmap -O -K >> $PS # auto correlation
awk '{print $3,$2,$4}' ~/Research/AuSREM-CM-moho/H-Mfiles/H-rtomo.xyz | gmt6 psxy -W.001 -Ss.18 -Cdevon_m.cpt $J $Rmap -O -K >> $PS # tomography
awk '{print $3,$2,$4}' ~/Research/AuSREM-CM-moho/H-Mfiles/*Recf*.xyz | gmt6 psxy -W.001 -Ss.18 -Cdevon_m.cpt $J $Rmap -O -K >> $PS # receiver functions

############### Active
awk '{print $3,$2,$4}' ~/Research/AuSREM-CM-moho/H-Mfiles/*reflect*.xyz | gmt6 psxy -W.001 -St.19 -Cdevon_m.cpt $J $Rmap -O -K >> $PS
awk '{print $3,$2,$4}' ~/Research/AuSREM-CM-moho/H-Mfiles/*refract*.xyz | gmt6 psxy -W.001 -St.19 -Cdevon_m.cpt $J $Rmap -O -K >> $PS


# echo 134 -30.5 Gawler Craton | gmt6 pstext $Rmap $J -F+f8.5p,Helvetica,firebrick -O -P -K  >> $PS
# echo 132 -26.5 Musgrave Block | gmt6 pstext $Rmap $J -F+f6.5p,Helvetica,darkcyan -O -P -K  >> $PS
# echo 141 -32.1 Curnamona Block | gmt6 pstext $Rmap $J -F+f5p,Helvetica,darkcyan -O -P -K  >> $PS
# echo 140.2 -28.7 Lake Eyre Basin | gmt6 pstext $Rmap $J -F+f6.5p,Helvetica,DarkGoldenrod -O -P -K  >> $PS
# echo 137.7 -30.5 THZ | gmt6 pstext $Rmap $J -F+f7.5p,Helvetica,DarkGreen+a-60 -O -P -K  >> $PS
# echo 138.6 -32 Flinders Ranges | gmt6 pstext $Rmap $J -F+f7.5p,Helvetica+a-80 -Gwhite@40 -O -P -K  >> $PS


gmt6 pscoast $Rmap $J -Bx -By -Na/.05p -A10 -P -Swhitesmoke  -Di -O -K -W.1p >> $PS #-A10+l -Ia

gmt6 psscale -Dx2c/-.5c+w3.2c/.25c+e+h -O -Cbil.cpt -K -Bx250 -By+l"Bouguer Gravity Anomaly (@~\225m@~m.s@+-2@+)" -P -V >> $PS
gmt6 psscale -Dx11.6c/3c+w4c/.3c+e -O -Cdevon_m.cpt -K -Bx5 -By+l"Moho (km)" -P -V >> $PS

# gmt6 psscale -Dx2c/-.7c+w3c/.3c+h -O -Cdevon.cpt -K -Bx.5 -By+l"Pms (sec)" -P -V >> $PS
###
##############
# gmt6 pscoast -R110E/155E/44S/9S -JM1.5/1i -P -O -Ba10f5 -Bwsne -Wfaint -N2/.05p  -EAU+gbisque -Gbrown -Sazure1 -Da -K -X10 -Y-5 -V --FORMAT_GEO_MAP=dddF >> $PS
# gmt6 psxy -R -J -O -T -X10 -Y-5 -V >> $PS
# AQT_1Q

gmt6 pslegend -Dx1.5c/.75c+w2c/1.3c+o-1c/-.5c -F+gwhite+p.1 -O $J $Rmap << EOF >> $PS
S 0.2c c 0.3c - black 0.2i 5G-6K
S 0.2c s 0.3c - black 0.2i Passive
S 0.2c t 0.3c - black 0.2i Active

EOF

####LEGEND
# gmt6 pslegend -DJLT+w2.7c/.4c+o-2.7c/-.4c -F+gwhite+p.1 -O $J $Rmap << EOF >> $PS
#S 0.2c s 0.3c gold - 0.2i 1Ds-1Dp
# EOF

#pstext $J $Rmap -P -O -K <<End >> $PS
#-40 -7 18 0 31 1 (b)
#End

# gmt6 psconvert -A -Tf -P -Z -Vq $PS
gmt6 ps2raster -A -E800 -Tj -P -Z -Vq $PS
open fig_grav.jpg
# cp fig1.pdf ~/Dropbox/NA_temp/moho_section/
