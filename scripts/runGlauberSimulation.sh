#!/bin/bash

#User-Defined Input Variables 
echo "Number of events: "
read N
#N=1000000
A=197
a=0.535
sigmaNN=33.1
option=makehists
optiona=progress


#Run the Glauber Model With the Above Options
root -l ../src/RunGlauberSimulation.cxx\($N\,$A\,$a\,$sigmaNN\,\"$option\"\,\"$optiona\"\)