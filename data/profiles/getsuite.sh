tartarusdir='/home/greg/Projects/LESbrary.jl/idealized/'
basename='three_layer_constant_fluxes_linear_hr'
hours=$1
suite="${hours}_hour_suite"

mkdir $suite
mkdir $suite/1m
mkdir $suite/2m
mkdir $suite/4m

scp -r tartarus:"${tartarusdir}${basename}${hours}*256*/*statistics.jld2" $suite/1m
scp -r tartarus:"${tartarusdir}${basename}${hours}*128*/*statistics.jld2" $suite/2m
scp -r tartarus:"${tartarusdir}${basename}${hours}*64*/*statistics.jld2"  $suite/4m

