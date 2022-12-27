#!/bin/bash
# Read spherical relocation (s3dloc.lis) file and plot the catalog
gmtset PAPER_MEDIA A3
gmtset LABEL_FONT_SIZE = 12p
gmtset BASEMAP_TYPE plain
gmtset ANNOT_FONT_SIZE_PRIMARY = 10p
gmtset ANNOT_FONT_SIZE_SECONDARY = 6p
gmtset HEADER_FONT_SIZE = 14p

rm *.ps

makecpt -Cseis -T-2/2/0.5 -I -Z > depth.cpt
makecpt -Cgray -T0/10000/1000 -I -Z > topo1.cpt
makecpt -Cwysiwyg -T-10/50/5 -I -Z > eqfocal.cpt


#
tit1="Time-distance diagram"
#tit2="Relocated"
# ======== Setting for out file ============
#out="compare_$if1.ps"
out="00.phases_compare_cut"
# Print setting
topodir=`pwd`

region=0/1200/0/200
# Plot medicator catalog

psbasemap -R$region -Glightyellow  -Ba200f2000g200:"Epicenter distance (km)":/2a50f300g50:"Traveltime (sec)":WeSn -JX6i/6i -K -V -X3i -X3i > $out.ps
gawk 'NR>=2&&$16 == 1 {print $13, $12}' init_data | psxy -R -JX -B -G200/0/0 -Sc0.1 -K -O -V >> $out.ps 
#gawk '{print $3, $2}' 00_Pn | psxy -R -JX -B -G0/200/0 -Sc0.1 -K -O -V >> $out.ps
gawk 'NR>=2&&$16 == 3 {print $13, $12}' init_data | psxy -R -JX -B -G0/0/200 -Sc0.1 -K -O -V >> $out.ps 
#gawk '{print $3, $2}' 00_Sn | psxy -R -JX -B -G0/0/0 -Sc0.1 -K -O -V >> $out.ps 
#   Number of EQ statistics within 10km
gawk 'NR>=2&&$16 == 1 {print $13, $12}' output_data_cut | psxy -R -JX -B -G0/0/0 -Sc0.1 -K -O -V >> $out.ps 
#gawk '{print $3, $2}' 00_Pn | psxy -R -JX -B -G0/200/0 -Sc0.1 -K -O -V >> $out.ps
gawk 'NR>=2&&$16 == 3 {print $13, $12}' output_data_cut | psxy -R -JX -B -G0/200/0 -Sc0.1 -K -O -V >> $out.ps 
#gawk '{print $3}' $if1 | pshistogram -Ba10f50:"Depth (km)":/1a100f500:"Numer of EQ":WSne \
# -R0/50/0/550 -JX4i/2i -Ggray -O -K -V  -L0.1p -Z0 -W10 -Y-6i -X6i >> $out.ps1

pslegend -R30/300/0/0.5 -Jx1i/-1i -O -K -V -Dx5.5/15/1.8i/1.6i/TC -Glightyellow -Y-0.2i -X-1i -Fthick >> $out.ps << END
N 7
G 0.5c
S 0.2i c 0.1i 200/0/0 1p,200/0/0 0.5i Pg phase input
G 0.5c
S 0.2i c 0.1i 0/0/200 1p,0/0/200 0.5i Pn phase input
G 0.5c
S 0.2i c 0.1i 0/0/0 1p,0/0/0 0.5i Pg phase output
G 0.5c
S 0.2i c 0.1i 0/200/0 1p,0/200/0 0.5i Pn phase output
END

# Clean up
rm *.cpt
#gs $out.ps
#ps2pdf $out.ps 
#ps2pdf $out.ps1
convert -trim -rotate 90 $out.ps $out.png
display $out.png
