#! /bin/sh
stamp=$(date +%s)
mkdir logs/${stamp}
cp *.* logs/${stamp}
