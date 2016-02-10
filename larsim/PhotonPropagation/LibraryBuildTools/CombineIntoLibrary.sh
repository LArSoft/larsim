#!/bin/bash

BatchID=$1
echo "Assembling library for $BatchID"

#Produce the list of files to combine
ls root/Files$BatchID/ > temp${BatchID}_FileList.txt

#Make a string to feed to roor
echo ".L AssembleSingleFile.C" >> rootstring.txt
echo "AssembleSingleFile(\"temp${BatchID}_FileList.txt\", \"root/Files$BatchID/\",\"photlibrary/lib$BatchID.root\")" >> rootstring.txt

#Make root execute this command
cat rootstring.txt | root -l

#clean up temporaryfiles
rm temp${BatchID}_FileList.txt
rm rootstring.txt
