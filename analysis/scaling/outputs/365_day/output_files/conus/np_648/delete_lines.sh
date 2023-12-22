#!/bin/bash

NEWOUTPUT="sm_parallel_conus_864.out"
## remove any line that contains Over 1 ##
sed -i '/Over 1/d' $NEWOUTPUT
## remove any line that contains in the Tsfc ##
sed -i '/in the Tsfc/d' $NEWOUTPUT
## remove any line that contains usually means ##
sed -i '/usually means/d' $NEWOUTPUT
## remove any line that contains atmopheric forcing ##
sed -i '/atmopheric forcing/d' $NEWOUTPUT
## remove any line that contains in snowmodel.par is set too low ##
sed -i '/in snowmodel.par is set too low/d' $NEWOUTPUT
## remove any line that contains Tsfc not converged ##
sed -i '/Tsfc not converged/d' $NEWOUTPUT
## remove any line that contains 1 m ##
sed -i '/1 m/d' $NEWOUTPUT
## remove any line that contains ZEROING OUT ##
sed -i '/ZEROING OUT/d' $NEWOUTPUT
## remove any line that contains ZEROING OUT ##
sed -i '/Used memory/d' $NEWOUTPUT
## remove any line that contains ZEROING OUT ##
sed -i '/Running: /d' $NEWOUTPUT