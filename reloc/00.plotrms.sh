#!/bin/bash
# Read spherical relocation (s3dloc.lis) file and plot the catalog
gmtset PAPER_MEDIA A3
gmtset LABEL_FONT_SIZE = 12p
gmtset BASEMAP_TYPE plain
gmtset ANNOT_FONT_SIZE_PRIMARY = 10p
gmtset ANNOT_FONT_SIZE_SECONDARY = 6p
gmtset HEADER_FONT_SIZE = 14p

#rm *.ps

#makecpt -Cseis -T1/5/0.5 -I -Z > num.cpt
#makecpt -Cgray -T0/10000/1000 -I -Z > topo1.cpt
#makecpt -Cwysiwyg -T-10/50/5 -I -Z > eqfocal.cpt

#echo -15   0   0 150 -5   0   0 150 >  cor.cpt
#echo  -5   0   0 150 -2   0   0 255 >> cor.cpt
#echo  -2   0   0 255  0 255 255 255 >> cor.cpt
#echo   0 255 255 255  2 255   0   0 >> cor.cpt
#echo   2 255   0   0  5 150   0   0 >> cor.cpt
#echo   5 150   0   0 15 150   0   0 >> cor.cpt
# setting for in file. Commend it when not use
out="00_plotrms"
# ======== RMS out ===================
psbasemap -R0.5/5.5/0.0/1.0 -Ba1f6g1:"Iteration number":/2a0.1f1.0g0.10:"RMS":WeSn -JX5i/5i  -K -V -X1i -Y3i > $out.ps
#
gawk '{print $1, $3}' s3dloc.rms | psxy -R -JX -B -Wthick,0/0/200,- -K -O -V >> $out.ps
gawk '{print $1, $3}' s3dloc.rms | psxy -R -JX -B -G0/0/200 -Sc0.08i -K -O -V >> $out.ps
#
#
#pslegend -R -JX -D6.2/2.0/1.5i/2.5i/TL -Fthick -K -O << END >> $out.ps
#N 1
#S 0.2c c 0.3c 200/0/200 0.25p 0.6c damp=0.01
#G 0.2c
#END
#
## =========== Number event reduce =============
psbasemap -R0.0/6.0/0/100 -Ba1f6g1:"Iteration number":/2a50f1500g100:"Number of event used each itteration":WeSn -JX5i/5i  -K -O -V -X8i -Y0i >> $out.ps
#
gawk '{print $1, $2}' s3dloc.rms | psxy -R -JX -B -Wthick,0/0/200,- -K -O -V >> $out.ps
gawk '{print $1, $2}' s3dloc.rms | psxy -R -JX -B -G0/0/200 -Sc0.08i -K -O -V >> $out.ps

#
#pslegend -R -JX -D6.2/1800/1.5i/2.5i/TL -Fthick -K -O << END >> $out.ps
#N 1
#S 0.2c c 0.3c 200/0/200 0.25p 0.6c damp=0.01
#G 0.2c
#END
#
# Clean up
rm *.cpt
#gs $out
#ps2pdf $out
convert -trim -rotate 90 $out.ps $out.png
display $out.png
