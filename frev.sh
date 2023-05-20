#! /bin/sh

dT=4
nsteps=40
div=100
N=3
r0=0.01
ic="circle"
flops=1

while getopts T:n:d:N:i:r: opt; do
    case $opt in
        T ) dT=${OPTARG};;
        n ) nsteps=${OPTARG};;
        d ) div=${OPTARG};;
        N ) N=${OPTARG};;
        i ) ic=${OPTARG};;
        r ) r0=${OPTARG};;
        x ) flops=${OPTARG};;
    esac
done

echo " "
echo "dT = $dT"
echo "nsteps = $nsteps"
echo "div = $div"
echo "N = $N"
echo "ic = $ic"
echo "r0 = $r0"
echo "flops = $flops"

echo " "
echo "init"
julia kinetic/init.jl -N $N -ic $ic
echo "sleep 1s"
sleep 1s

echo "sim"
julia kinetic/kinetic.jl -dT $dT -nsteps $nsteps -div $div -r0 $r0
echo "rev"
julia kinetic/frev.jl
echo "sim"
julia kinetic/kinetic.jl -dT $dT -nsteps $nsteps -div $div -r0 $r0

echo "animate"
python3 kinetic/animate.py
echo "view"
viewnior out.gif
