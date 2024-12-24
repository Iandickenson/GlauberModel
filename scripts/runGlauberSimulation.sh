#!/bin/bash

#Input Variables 
N=100
A=197
a=0.535
sigmaNN=33.1
option=makehists


#Run the Glauber Model With the Above Options
root -l ../src/RunGlauberSimulation.cxx\($N\,$A\,$a\,$sigmaNN\,\"$option\"\)