#!/bin/bash

#Real job
jobsub --opportunistic --X509_USER_PROXY /scratch/bjpjones/grid/bjpjones.uboone.proxy -g -N 9375 -dOUT /uboone/data/users/bjpjones/OpticalProduction  -q OpticalLibraryBuild_Grid.sh `whoami` `pwd`
#jobsub [Client Options] user_script [user_script_args]
#--opportunistic
#--X509_USER_PROXY=/scratch/bjpjones/grid/bjpjones.uboone.proxy
#-g = grid
#-N = number of jobs
#-dOUT dir = output directory
#-q = mail_on_error
#user_script = OpticalLibraryBuild_Grid.sh
#`whoami`=user_name
#`pwd`=working directory

#Partial job
#jobsub --opportunistic --X509_USER_PROXY /scratch/bjpjones/grid/bjpjones.uboone.proxy -g -N 1500 -dOUT /uboone/data/users/bjpjones/OpticalProduction  -q OpticalLibraryBuild_Grid.sh `whoami` `pwd`

#Test job
#jobsub --opportunistic --X509_USER_PROXY /scratch/bjpjones/grid/bjpjones.uboone.proxy -g -N 5 -dOUT /uboone/data/users/bjpjones/OpticalProduction  -q OpticalLibraryBuild_Grid.sh `whoami` `pwd`
