#! /bin/sh
stamp=$(date +%s)
mkdir logs/${stamp}
cp *.* logs/${stamp}
mkdir logs/${stamp}/kinetic
cp kinetic/*.* logs/${stamp}/kinetic
