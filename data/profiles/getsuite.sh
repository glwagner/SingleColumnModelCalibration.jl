#!/usr/bin/env bash

shopt -s extglob

tartarus='greg@tartarus.mit.edu:/home/greg/Projects/LESbrary.jl/idealized/'
basename='three_layer_constant_fluxes_linear_hr'

hours=6
suite="${hours}_hour_suite"

mkdir $suite
mkdir $suite/1m
mkdir $suite/2m
mkdir $suite/4m

scp -r ${tartarus}${basename}${hours}*256*/*statistics.jld2 $suite/1m/ 
scp -r ${tartarus}${basename}${hours}*128*/*statistics.jld2 $suite/2m/
scp -r ${tartarus}${basename}${hours}*64*/*statistics.jld2  $suite/4m/

hours=12
suite="${hours}_hour_suite"

mkdir $suite
mkdir $suite/1m
mkdir $suite/2m
mkdir $suite/4m

scp -r ${tartarus}${basename}${hours}*256*/*statistics.jld2 $suite/1m/ 
scp -r ${tartarus}${basename}${hours}*128*/*statistics.jld2 $suite/2m/
scp -r ${tartarus}${basename}${hours}*64*/*statistics.jld2  $suite/4m/

hours=24
suite="${hours}_hour_suite"

mkdir $suite
mkdir $suite/1m
mkdir $suite/2m
mkdir $suite/4m

scp -r ${tartarus}${basename}${hours}*256*/*statistics.jld2 $suite/1m/ 
scp -r ${tartarus}${basename}${hours}*128*/*statistics.jld2 $suite/2m/
scp -r ${tartarus}${basename}${hours}*64*/*statistics.jld2  $suite/4m/

hours=48
suite="${hours}_hour_suite"

mkdir $suite
mkdir $suite/1m
mkdir $suite/2m
mkdir $suite/4m

scp -r ${tartarus}${basename}${hours}*256*/*statistics.jld2 $suite/1m/ 
scp -r ${tartarus}${basename}${hours}*128*/*statistics.jld2 $suite/2m/
scp -r ${tartarus}${basename}${hours}*64*/*statistics.jld2  $suite/4m/

hours=72
suite="${hours}_hour_suite"

mkdir $suite
mkdir $suite/1m
mkdir $suite/2m
mkdir $suite/4m

scp -r ${tartarus}${basename}${hours}*256*/*statistics.jld2 $suite/1m/ 
scp -r ${tartarus}${basename}${hours}*128*/*statistics.jld2 $suite/2m/
scp -r ${tartarus}${basename}${hours}*64*/*statistics.jld2  $suite/4m/
