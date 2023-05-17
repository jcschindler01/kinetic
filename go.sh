#! /bin/sh

# usage: go.sh -N 10 -T dT -n nsteps -d div 

dT=2
nsteps=40
div=100
N=100
r0=0.01
ic="chain"

while getopts T:n:d:N:i:r: opt; do
    case $opt in
        T ) dT=${OPTARG};;
        n ) nsteps=${OPTARG};;
        d ) div=${OPTARG};;
        N ) N=${OPTARG};;
        i ) ic=${OPTARG};;
        r ) r0=${OPTARG};;
    esac
done

echo " "
echo "dT = $dT"
echo "nsteps = $nsteps"
echo "div = $div"
echo "N = $N"
echo "ic = $ic"
echo "r0 = $r0"

echo " "
echo "init"
julia kinetic/init.jl -N $N -ic $ic
echo "sleep 1s"
sleep 1s
echo "sim"
julia kinetic/kinetic.jl -dT $dT -nsteps $nsteps -div $div -r0 $r0
echo "animate"
python3 kinetic/animate.py
echo "view"
viewnior out.gif
