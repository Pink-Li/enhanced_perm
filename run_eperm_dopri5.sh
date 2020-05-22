#!/bin/bash
clear
make -f make_dopri5 clean --silent
make -f make_dopri5 > make.log
ulimit -s unlimited 
./eperm_dopri5
#cd src
#matlab -nosplash -nodesktop -logfile 'matlab.log' -r result_viewer_1169
#cd ..
#reset


