#!/bin/bash
clear
make -f make_odepack clean --silent
make -f make_odepack > make.log
ulimit -s unlimited 
./eperm_odepack
#cd src
#matlab -nosplash -nodesktop -logfile 'matlab.log' -r result_viewer_1169
#cd ..
#reset


