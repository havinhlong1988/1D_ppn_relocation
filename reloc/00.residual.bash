#!/bin/bash
# PARAMETER SETTING -------------------------------
gmtset ANNOT_FONT_SIZE 14p 
gmtset ANNOT_OFFSET 0.05 
gmtset LABEL_FONT_SIZE 16p 
gmtset LABEL_OFFSET 0.03 
gmtset FRAME_PEN 1.5p 
gmtset FRAME_WIDTH 0.06i 
gmtset TICK_LENGTH 3.0p 
gmtset TICK_PEN 1.5p


fn=00_residual

# DATA PROCESSING ---------------------------------
gawk '{print $2}' ../1.ppn-it5-/residual_001 > 1dinput.plot
gawk '{print $2}' ../1.ppn-it5-/residual_005 > 1doutput.plot
#gawk '$2 >=-5.0 || $2 <= 5.0 {print $2}' residual_001 > ite1.plot
#gawk '$2 >=-5.0 || $2 <= 5.0 {print $2}' residual_002 > ite2.plot
#gawk '$2 >=-5.0 || $2 <= 5.0 {print $2}' residual_003 > ite3.plot
gawk '{print $2}' residual_001 > ite1.plot
gawk '{print $2}' residual_002 > ite2.plot
gawk '{print $2}' residual_003 > ite3.plot
gawk '{print $2}' residual_004 > ite4.plot
gawk '{print $2}' residual_005 > ite5.plot

#gawk "NR>=3 && $1==121.32 && $2==23.94 {print $4,$3}" vpvstomo.dat > 1dvp1.plot
#gawk "NR>=3 && $1==121.32 && $2==23.94 {print $7,$3}" vpvstomo.dat > 1dvs1.plot
#gawk "NR>=3 && $1==121.32 && $2==23.94 {print $10,$3}" vpvstomo.dat > 1dvr1.plot



# PLOTTING ----------------------------------------
#1d vp & vs models
psbasemap -Glightyellow -B1a1f50g1:"Residual (sec)":/a100g50f12000:"Number of reading":WeSn -JX10i/10i -R-3.0/3.0/0/80 -K -V -Y2.0 -X2.0 > $fn.ps
#sxy zome-in.plot -R -JX -B -W2/0 -G235 -K -O -V >> $fn.ps
#psxy 1dmods_vp.dat -R -JX -W2/150 -MX -K -O -V >> $fn.ps
#pshistogram 1dinput.plot -Ba:"Depth (km)":/100:"Numer of EQ":WSne -R-3.0/3.0/0/2000 -JX6i/6i -Gred -K -V -O -L0.1p -Z0 -W1 >> $fn.ps
pshistogram 1dinput.plot -R -JX -Z0 -W0.05 -S -L2.0p,red,. -K -O -V >> $fn.ps
pshistogram 1doutput.plot -R -JX -Z0 -W0.05 -S -L2.0p,blue,. -K -O -V >> $fn.ps
pshistogram ite1.plot -R -JX -Z0 -W0.05 -S -L1.0p,red -K -O -V >> $fn.ps
pshistogram ite2.plot -R -JX -Z0 -W0.05 -S -L1.0p,green -K -O -V >> $fn.ps
pshistogram ite3.plot -R -JX -Z0 -W0.05 -S -L1.0p,blue -K -O -V >> $fn.ps
pshistogram ite4.plot -R -JX -Z0 -W0.05 -S -L1.0p,black -K -O -V >> $fn.ps
pshistogram ite5.plot -R -JX -Z0 -W0.05 -S -L1.0p,gray -K -O -V >> $fn.ps
#
pslegend -R30/300/0/0.5 -Jx1i/-1i -O -K -V -Dx5.5/25/1.8i/3.0i/TC -Glightyellow -Y-0.2i -X-1i -Fthick >> $fn.ps << END
N 7
G 0.5c
S 0.1i c 0.1i 200/0/0 1p,200/0/0 0.5i input 1D
G 0.5c
S 0.1i c 0.1i 0/0/200 1p,0/0/200 0.5i output 1D
G 0.5-
S 0.2i - 0.3i 200/0/0 1p,200/0/0 0.5i reloc it1
G 0.5-
S 0.2i - 0.3i 0/200/0 1p,0/200/0 0.5i reloc it2
G 0.5-
S 0.2i - 0.3i 0/0/200 1p,0/0/200 0.5i reloc it3
G 0.5-
S 0.2i - 0.3i 0/0/0 1p,0/0/0 0.5i reloc it4
G 0.5-
S 0.2i - 0.3i 128/128/128 1p,128/128/128 0.5i reloc it5
END


# CLEANING ---------------------------------------------
#rm -f *.plot

# DISPLAYING -------------------------------------------
convert -trim -rotate 90 $fn.ps $fn.png
display $fn.png

