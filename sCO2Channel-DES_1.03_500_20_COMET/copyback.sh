#!/bin/bash
cd $1
taskid=$3
tar -cf $2/node$taskid.tar processor*
