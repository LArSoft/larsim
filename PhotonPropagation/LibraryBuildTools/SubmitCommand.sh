#!/bin/bash

#Real job
jobsub --opportunistic --X509_USER_PROXY /scratch/bjpjones/grid/bjpjones.uboone.proxy -g -N 9375 -dOUT /uboone/data/users/bjpjones/OpticalProduction  -q OpticalLibraryBuild_Grid.sh `whoami` `pwd`

#Partial job
#jobsub --opportunistic --X509_USER_PROXY /scratch/bjpjones/grid/bjpjones.uboone.proxy -g -N 1500 -dOUT /uboone/data/users/bjpjones/OpticalProduction  -q OpticalLibraryBuild_Grid.sh `whoami` `pwd`

#Test job
#jobsub --opportunistic --X509_USER_PROXY /scratch/bjpjones/grid/bjpjones.uboone.proxy -g -N 5 -dOUT /uboone/data/users/bjpjones/OpticalProduction  -q OpticalLibraryBuild_Grid.sh `whoami` `pwd`
