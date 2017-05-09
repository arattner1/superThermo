#!/bin/bash
NNODES=$1
NTASKPERNODE=$2

# Make directory for initial backup
mkdir -p BACK0

# Tar up the common directories
FILEN=common.tar
if [ -f $FILEN ]
then
 echo "$FILEN already exists"
else
 tar -cvf $FILEN 0 system constant
fi

pnumber=0
for (( nn=1; nn<=$NNODES; nn++ ))
 do 
  FILEN=node$nn.tar
  if [ -f $FILEN ]
  then
    echo "$FILEN already exists"
  else 
    mkdir -p TPROC 
    cd TPROC
    for (( mm=1; mm<=$NTASKPERNODE; mm++ ))
     do
      mv ../processor$pnumber .
      ((pnumber+=1))
     done
    tar -cvf $FILEN processor*
    mv $FILEN ../
    mv processor* ../BACK0
    cd ../
  fi
 done
  
