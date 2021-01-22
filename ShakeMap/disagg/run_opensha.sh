#!/bin/bash

set -o errexit

if [[ $# -lt 1 ]];then
	echo "must supply class to run (and optional arguments)"
	exit 2
fi

java=`which java`

classpath="opensha-cybershake-all.jar"
$java -Xms512M -Xmx2G -cp $classpath $@
