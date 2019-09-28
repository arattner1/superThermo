#!/bin/bash
cd $1
tar -xf $2/common.tar
taskid=$3
tar -xf $2/node$taskid.tar
