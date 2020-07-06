#!/bin/bash
if [ $# != 2 ] ; then
	echo "cat Prog_1809.py |./ffp.sh 255 308"
	exit
fi
from=$(($1))
to=$(($2))
tail -n +$from |head -n $(($to-$from))
