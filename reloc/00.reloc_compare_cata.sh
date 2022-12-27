#!/bin/bash
# Read spherical relocation (s3dloc.lis) file and plot the catalog
gmtset PAPER_MEDIA A3
gmtset LABEL_FONT_SIZE = 12p
gmtset BASEMAP_TYPE plain
gmtset ANNOT_FONT_SIZE_PRIMARY = 10p
gmtset ANNOT_FONT_SIZE_SECONDARY = 6p
gmtset HEADER_FONT_SIZE = 14p

rm *.ps
makecpt -Cseis -T0/50/10 -I -Z > eqfocal.cpt
makecpt -Cseis -T-2/2/0.5 -I -Z > depth.cpt
makecpt -Cgray -T0/10000/1000 -I -Z > topo1.cpt
#makecpt -Cwysiwyg -T-10/50/5 -I -Z > eqfocal.cpt

# Input file
if="00_evt"
if1="00_evt_cut.out"
if2="fault_lv1.txt"
if3="00_sta"
rayfile="ray_path_003"
#   ===========  Check file length =============== 
(wc -l $if) | gawk '{print $1}' > temp 
fl=`awk '{print $1}' temp`
#echo $fl
(wc -l $if1) | gawk '{print $1}' > temp1
fl1=`awk '{print $1}' temp1`
#echo $fl1
#   =========== Check the moho depth ===============
dm=`gawk 'FNR==5{print $1}' input`
echo $dm
#
tit1="Origin catalog"
tit2="Relocated"
# ======== Setting for out file ============
out="00_compare_$if1"
out1="00_combine_$if1"
out2="00_pg_$if1" # Pg output file
out3="00_pn_$if1" # Pg output file
# Print setting
topodir=`pwd`
#   =====================================================================================================
#   Produce the raypath file to plot by fortran program producerayfile.f90
gfortran -o producerayfile producerayfile.f90
echo $rayfile > tmp.in
./producerayfile < tmp.in 
#$rayfile
#   =====================================================================================================
region=100.0/110.0/17.5/25.5

# ================= Input catalog ============================
grdgradient $topodir/topo.grd -G$topodir/topo.shd -A350 -Ne1.0 -M
grdimage -JM5i $topodir/topo.grd -R$region -I$topodir/topo.shd  -C$topodir/topo1.cpt \
-Ba1f1:."$if-$tit1":WeNs -K -V -Y15.0 -X2.5 > $out.ps
#
psbasemap -J -R -B -P -K -O -V >> $out.ps
pscoast -J -R -B -Df -Na/4 -W1p/0 -K -O -V >> $out.ps
# fault
gawk '{print $1,$2}' $if2 | psxy -J -R -MX -W2.5/red -O -K -V  >> $out.ps
# station
gawk '{print $3,$2,0.0,$1}' $if3 | psxyz -R -J -B -W1/0,black -G0/0/255 -E180/90 -St0.1i -K -O -V >> $out.ps
# EQ focal
gawk '{print $4,$3,$5, $6*0.02}' $if | psxy -J -R -B -Ceqfocal.cpt -W1/0,black -Sci -Wthin -K -O >> $out.ps
psscale -Ceqfocal.cpt -D9.4/1.5/-1.0i/0.05i -L -K -O -V >> $out.ps
echo 107.0 19.2 12 90 3 Mc 'Depth (km)' | pstext -J -R -G0 -S8/255 -E180/90 -K -O -V >> $out.ps
#  ---------------- x cross -------------------------
psbasemap -R100/110/-5/50 -JX5i/-1.5i -B1:"Long (deg)":/1a10f50g10WnSe:"Depth (km)": -P -K -O -V -X0i -Y-1.7i >> $out.ps
gawk '{print $3,-$4/1000}' $if3 | psxy -J -R -B  -St0.05i -W1/0,blue -G0/0/200 -K -O >> $out.ps
gawk '{print $4,$5,$5,$6*0.02}' $if | psxy -J -R -B -Ceqfocal.cpt -W2/0,black -Sci -K -O -V >> $out.ps
# Moho line
psxy -: -R -JX -MX -W4/,blue,- -K -O -V >> $out.ps << END
$dm 99
$dm 111
END
echo 102 40 12 0 3 ML Moho=$dm km | pstext -JX -R -Gblue -S5/255 -N -K -O -V >> $out.ps
#  ----------------- y cross ------------------------------
psbasemap -R-5/50/18/25.5  -JX1.5i/4.25i -B1a10f50g10wNsE:"Depth (km)":/1:"Lat (deg)": -P -K -O -X5.2i -Y1.7i >> $out.ps
gawk '{print -$4/1000, $2}' $if3 | psxy -J -R -B  -St0.05i -W1/0,blue -G0/0/200 -K -O >> $out.ps
gawk '{print $5,$3,$5,$6*0.02}' $if | psxy -J -R -B -Ceqfocal.cpt -W2/0,black -Sci -K -O -V >> $out.ps
# Moho line
psxy -: -R -JX -MX -W4/,blue,- -K -O -V >> $out.ps << END
10 $dm
30 $dm
END
echo 40 18 12 90 3 ML Moho=$dm km | pstext -JX -R -Gblue -S5/255 -N -K -O -V >> $out.ps


# ============== Relocated catalog ==========================
grdimage -JM5i $topodir/topo.grd -R$region -I$topodir/topo.shd  -C$topodir/topo1.cpt \
-Ba1f1:."$if1-$tit2":WNse -K -O -V -Y0.0 -X7.0 >> $out.ps
psbasemap -J -R -B -P -K -O -V >> $out.ps
pscoast -J -R -B -Df -Na/4 -W1p/0 -K -O >> $out.ps
# fault
gawk '{print $1,$2}' $if2 | psxy -J -R -MX -O -K -V -W2.5/red >> $out.ps
# station
gawk '{print $3, $2}' $if3 | psxy -R -JM -W4/0,black -St0.1i -Gblue -K -O -V >> $out.ps
# EQ focal
gawk '{print $4,$3,$5, $6*0.02}' $if1 | psxy -J -R -B -Ceqfocal.cpt -W1/0,red -Sci -Wthin -K -O >> $out.ps
psscale -Ceqfocal.cpt -D9.4/1.5/-1.0i/0.05i -L -K -O -V >> $out.ps
echo 107.0 19.2 12 90 3 Mc 'Depth (km)' | pstext -J -R -G0 -S8/255 -E180/90 -K -O -V >> $out.ps
#  ---------------- x cross -------------------------
psbasemap -R100/110/-5/50 -JX5i/-1.5i -B1:"Long (deg)":/1a10f50g10WnSe:"Depth (km)": -P -K -O -V -X0i -Y-1.7i >> $out.ps
gawk '{print $3,-$4/1000}' $if3 | psxy -J -R -B  -St0.05i -W1/0,blue -G0/0/200 -K -O >> $out.ps
gawk '{print $4,$5,$5,$6*0.02}' $if1 | psxy -J -R -B -Ceqfocal.cpt -W2/0,red -Sci -K -O -V >> $out.ps
# Moho line
psxy -: -R -JX -MX -W4/,blue,- -K -O -V >> $out.ps << END
$dm 100
$dm 110
END
echo 102 40 12 0 3 ML Moho=$dm km | pstext -JX -R -Gblue -S5/255 -N -K -O -V >> $out.ps
#  ----------------- y cross ------------------------------
psbasemap -R-5/50/18/25.5  -JX1.5i/4.25i -B2a10f50g10wNsE:"Depth (km)":/1:"Lat (deg)": -P -K -O -X5.2i -Y1.7i >> $out.ps
gawk '{print -$4/1000, $2}' $if3 | psxy -J -R -B  -St0.05i -W1/0,blue -G0/0/200 -K -O >> $out.ps
gawk '{print $5,$3,$5,$6*0.02}' $if1 | psxy -J -R -B -Ceqfocal.cpt -W2/0,red -Sci -K -O -V >> $out.ps
# Moho line
psxy -: -R -JX -MX -W4/,blue,- -K -O -V >> $out.ps << END
10 $dm
30 $dm
END
echo 40 18 12 90 3 ML Moho=$dm km | pstext -JX -R -Gblue -S5/255 -N -K -O -V >> $out.ps
# Legend 2
pslegend -R0/9/0/0.5 -Jx1i/-1i -O -K -V -Dx4.5/0/1.5i/1.5i/TC -Y-0.2i -X-1i -Fthick >> $out.ps << END
N 7
G 0.5c
S 0.2i c 0.2i white 1p,red 0.5i $tit2
G 0.5c
S 0.0i . 0.0i red 0.5p 0.5i $fl1 events
END

# Legend 1
pslegend -R0/9/0/0.5 -Jx1i/-1i -O -K -V -Dx4.5/0/1.5i/1.5i/TC -Y0.0i -X-7.85i -Fthick >> $out.ps << END
N 7
G 0.5c
S 0.2i c 0.2i white 1p,back 0.5i $tit1
G 0.5c
S 0.2i . 0.2i white 0.5p,red 0.5i $fl events
END
## Statistical grap
gawk '{print $5}' $if | pshistogram -Ba5f50:"Depth (km)":/1a100f500:"Numer of EQ"::."Input catalog histogram":WSne \
 -R0/40/0/450 -JX5i/2i -Gblack -K -V -O -L0.2p -Z0 -W1 -X-4.0i -Y-4.5i >> $out.ps

gawk '{print $5}' $if1 | pshistogram -Ba5f50:"Depth (km)":/1a50f500:"Numer of EQ"::."Relocated catalog histogram":WSne \
 -R0/40/0/150 -JX5i/2i -Gred -K -V -O -L0.2p -Z0 -W1 -X8.0i >> $out.ps


## ==============================   Combined plot ===============================================================

grdimage -JM7i $topodir/topo.grd -R$region -I$topodir/topo.shd  -C$topodir/topo1.cpt \
-Ba1f1:."$if1":WeNs -K -V -Y8.0 -X2.5 > $out1.ps
psbasemap -JM7i -R$region -B2/2wesn -P -K -O -V >> $out1.ps
pscoast -J -R -B -Df -Na/4 -I1/3/0/0/255 -W1p/0 -K -O >> $out1.ps
# station
gawk '{print $3,$2,0.0,$1}' $if2 | psxyz -R -J -B -W1/0,black -G0/0/255 -E180/90 -St0.15i -K -O -V >> $out1.ps
# catalog 1
gawk '{print $4,$3,$5, $6*0.025}' $if | psxy -J -R -B -Ceqfocal.cpt -W1/0,black -Sci -K -O >> $out1.ps
# catalog 2
gawk '{print $4,$3,$5, $6*0.025}' $if1 | psxy -J -R -B -Ceqfocal.cpt -W1/0,red -Sci -K -O >> $out1.ps
# Depth bar
psscale -Ceqfocal.cpt -D13/2.2/-1.5i/0.07i -L -K -O -V >> $out1.ps
echo 107.0 19.2 12 90 3 Mc 'Depth (km)' | pstext -J -R -G0 -S8/255 -E180/90 -K -O -V >> $out1.ps

#   x cross
psbasemap -R100/110/-5/50 -JX7i/-2i -B1:"Long (deg)":/2a10f100g10WnSe:"Depth (km)": -P -K -O -V -X0i -Y-2.2i >> $out1.ps
gawk '{print $3,-$4/1000}' $if3 | psxy -J -R -B  -St0.05i -W1/0,blue -G0/0/200 -K -O >> $out1.ps
gawk '{print $4,$5,$6*0.02}' $if | psxy -J -R -B -Sci  -W1/0,black -K -O -V >> $out1.ps
gawk '{print $4,$5,$6*0.02}' $if1 | psxy -J -R -B -Sci -W1/0,red -K -O -V >> $out1.ps
# Moho line
psxy -: -R -JX -MX -W4/,blue,- -K -O -V >> $out1.ps << END
$dm 100
$dm 110
END
echo 102 40 12 0 3 ML Moho=$dm km | pstext -JX -R -Gblue -S5/255 -N -K -O -V >> $out1.ps
#   y cross
psbasemap -R5/50/18/25.5  -JX2i/5.60i -B2a10f100g10wNsE:"Depth (km)":/1:"Lat (deg)": -P -K -O -X7.2i -Y2.2i >> $out1.ps
gawk '{print $5 ,$3 ,$6*0.025}' $if | psxy -J -R -B -Sci -W1/0,black -K -O >> $out1.ps
gawk '{print $5 ,$3 ,$6*0.025}' $if1 | psxy -J -R -B -Sci -W1/0,red -K -O >> $out1.ps
# Moho line
psxy -: -R -JX -MX -W4/,blue,- -K -O -V >> $out1.ps << END
10 $dm
30 $dm
END
echo 40 18 12 90 3 ML Moho=$dm km | pstext -JX -R -Gblue -S5/255 -N -K -O -V >> $out1.ps

# Legend
pslegend -R0/9/0/0.5 -Jx2i/-2i -O -Dx4.5/0/1.5i/1.5i/TC -Y-0.2i -X-1i -Fthick >> $out1.ps << END
N 5
G 0.5c
S 0.2i c 0.2i black 0.5p 0.3i $tit1
G 0.5c
S 0.2i . 0.2i black 0.5p 0.3i $fl events
G 0.5c
N 7
S 0.2i c 0.2i red 0.5p 0.3i $tit2
G 0.5c
S 0.2i . 0.2i red 0.5p 0.3i $fl1 events
END
#S 0.15i c 0.2i red 0.25p 0.3i $tit2
#	==========================================================================================================
#		Ray coverage 
#
#	Vertical ray coverage 
#awk -F " " '{gsub(/\#/,">",$1);print $1, $2, $3, $4}' $rayfile > raypath

#	==========================================================================================================
echo "--->  plot ray coverage Pg"
# =============== All ray coverage ====================
# Horizontal and horizontal ray
psbasemap -JM7i -R100/110/17.5/25.5 -B2a2f10:."Horizontaly ray coverage":/2a2f10WesN -K -V -X2i -Y5i > $out2.ps
pscoast -J -R -B -Df -Na/4 -W1p/0 -K -O -V >> $out2.ps
psxy pghray.plot -MX -R -J -W0.001,gray -K -O -V >> $out2.ps
gawk '{print $3, $2}' $if3 | psxy -R -JM -B -G0/0/200 -St0.05i -K -O -V >> $out2.ps
gawk '{print $4, $3}' $if1 | psxy -R -JM -B -G200/0/0 -Sc0.01i -K -O -V >> $out2.ps
#
# Vertical and horizontal ray
psbasemap -JX7i/-2i -R100/110/-5/50 -B2a2f10:"Long (deg)":/2a50f10g10:"Depth (km)":WeSn -P -K -O -V -Xi -Y-2.2i >> $out2.ps
psxy pgxray.plot -MX -R -J -W0.01,gray -K -O -V >> $out2.ps
#gawk ' $4==2 {print $2,$3}' raypath | psxy -m -R -J -B -W0.01,orange -K -O -V >> $out2.ps
gawk '{print $3, $4/-1000}' $if3 | psxy -m -R -J -B -G0/0/200 -St0.05i -K -O -V >> $out2.ps
gawk '{print $4, $5}' $if1 | psxy -m -R -J -B -G200/0/0 -Sc0.005i -K -O -V >> $out2.ps
psxy -: -R -JX -MX -W4/,blue,- -K -O -V >> $out2.ps << END
$dm 99
$dm 111
END
echo 102 40 12 0 3 ML Moho=$dm km | pstext -JX -R -Gblue -S5/255 -N -K -O -V >> $out2.ps
#
# Vertical and horizontal ray
psbasemap -JX2i/6i -R-5/50/18/25.5 -B2a10f50g10:."Depth (km)":/2a2f10g2:"Lat (deg)":wEsN -P -K -O -V -X7.2i -Y2.2i >> $out2.ps
psxy pgyray.plot -MX -R -J -W0.01,gray -K -O -V >> $out2.ps
#gawk ' $4==2 {print $3,$1}' raypath | psxy -m -R -J -B -W0.01,orange -K -O -V >> $out2.ps
gawk '{print $4/-1000, $2}' $if3 | psxy -m -R -J -B -G0/0/200 -St0.05i -K -O -V >> $out2.ps
gawk '{print $5, $3}' $if1 | psxy -m -R -J -B -G200/0/0 -Sc0.005i -K -O -V >> $out2.ps
psxy -: -R -JX -MX -W4/,blue,- -K -O -V >> $out2.ps << END
10 $dm
30 $dm
END
echo 40 18 12 90 3 ML Moho=$dm km | pstext -JX -R -Gblue -S5/255 -N -K -O -V >> $out2.ps

echo "--->  plot ray coverage Pn"
# =============== All ray coverage ====================
# Horizontal and horizontal ray
psbasemap -JM7i -R100/110/17.5/25.5 -B2a2f10:."Horizontaly ray coverage":/2a2f10WesN -K -V -X2i -Y3i > $out3.ps
pscoast -J -R -B -Df -Na/4 -W1p/0 -K -O -V >> $out3.ps
psxy pnhray.plot -MX -R -J -W0.001,orange -K -O -V >> $out3.ps
gawk '{print $3, $2}' $if3 | psxy -R -JM -B -G0/0/200 -St0.05i -K -O -V >> $out3.ps
gawk '{print $4, $3}' $if1 | psxy -R -JM -B -G200/0/0 -Sc0.01i -K -O -V >> $out3.ps
#
# Vertical and horizontal ray
psbasemap -JX7i/-2i -R100/110/-5/50 -B2a2f10:"Long (deg)":/2a50f10g10:"Depth (km)":WeSn -P -K -O -V -Xi -Y-2.2i >> $out3.ps
#gawk ' $4==1 {print $2,$1}' raypath | psxy -m -R -J -B -W0.01,gray -K -O -V >> $out3.ps
psxy pnxray.plot -MX -R -J -W0.01,orange -K -O -V >> $out3.ps
gawk '{print $3, $4/-1000}' $if3 | psxy -m -R -J -B -G0/0/200 -St0.05i -K -O -V >> $out3.ps
gawk '{print $4, $5}' $if1 | psxy -m -R -J -B -G200/0/0 -Sc0.005i -K -O -V >> $out3.ps
psxy -: -R -JX -MX -W4/,blue,- -K -O -V >> $out3.ps << END
$dm 99
$dm 111
END
echo 102 40 12 0 3 ML Moho=$dm km | pstext -JX -R -Gblue -S5/255 -N -K -O -V >> $out3.ps
#
# Vertical and horizontal ray
psbasemap -JX2i/6i -R-5/50/18/25.5 -B2a10f50g10:."Depth (km)":/2a2f10g2:"Lat (deg)":wEsN -P -K -O -V -X7.2i -Y2.2i >> $out3.ps
#gawk ' $4==1 {print $3,$1}' $rayfile | psxy -m -R -J -B -W0.01,gray -K -O -V >> $out3.ps
psxy pnyray.plot -MX -R -J -W0.01,orange -K -O -V >> $out3.ps
gawk '{print $4/-1000, $2}' $if3 | psxy -m -R -J -B -G0/0/200 -St0.05i -K -O -V >> $out3.ps
gawk '{print $5, $3}' $if1 | psxy -m -R -J -B -G200/0/0 -Sc0.005i -K -O -V >> $out3.ps
psxy -: -R -JX -MX -W4/,blue,- -K -O -V >> $out3.ps << END
10 $dm
30 $dm
END
echo 40 18 12 90 3 ML Moho=$dm km | pstext -JX -R -Gblue -S5/255 -N -K -O -V >> $out3.ps
# Clean up
rm depth.cpt eqfocal.cpt topo1.cpt *.plot
#gs $out
convert -trim -rotate 90 $out.ps $out.png 
convert -trim -rotate 90 $out1.ps $out1.png
convert -trim -rotate 90 $out2.ps $out2.png
convert -trim -rotate 90 $out3.ps $out3.png
#
display $out.png