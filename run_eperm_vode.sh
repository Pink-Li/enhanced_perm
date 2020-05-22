#!/bin/bash
clear
make -f make_vode clean --silent
make -f make_vode > make.log
ulimit -s unlimited 
./eperm_vode
#cd src
#matlab -nosplash -nodesktop -logfile 'matlab.log' -r result_viewer_1169
#cd ..
#reset


