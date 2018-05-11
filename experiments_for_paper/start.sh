#! /bin/bash

setupvals=(1 2)
nvals=(500 1000)
pvals=(6 12)
sigmavals=(0.1 1.0 2.0 5.0 10.0)
kvals=(2 3 4 5)
etavals=(0.1 0.2)
algvals=('R' 'X' 'S' 'T' 'U' 'RS')

# setup = as.numeric(args[1])
# n = as.numeric(args[2])
# p = as.numeric(args[3])
# sigma = as.numeric(args[4])
# k = as.numeric(args[5])
# eta = as.numeric(args[6])
# alg = as.numeric(args[7])
# NREP = as.numeric(args[8])

reps=200

for ((i1=0; i1<${#setupvals[@]} ;i1++))
do
for ((i2=0; i2<${#nvals[@]} ;i2++))
do
for ((i3=0; i3<${#pvals[@]} ;i3++))
do
for ((i4=0; i4<${#sigmavals[@]} ;i4++))
do
for ((i5=0; i5<${#kvals[@]} ;i5++))
do
for ((i6=0; i6<${#etavals[@]} ;i6++))
do
for ((i7=0; i7<${#algvals[@]} ;i7++))
do
    setup=${setupvals[$i1]}
    n=${nvals[$i2]}
    p=${pvals[$i3]}
    sigma=${sigmavals[$i4]}
    k=${kvals[$i5]}
    eta=${etavals[$i6]}
    alg=${algvals[$i7]}

    fnm="logging/progress-$setup-$n-$p-$sigma-$k-$eta-$alg-$reps.out"
    #echo $fnm

    Rscript run_simu.R $setup $n $p $sigma $k $eta $alg $reps 2>&1 | tee $fnm &
    #echo "Rscript run_simu.R $setup $n $p $sigma $k $eta $alg $reps 2>&1 | tee $fnm"
    #R CMD BATCH --no-save --no-restore "--args $setup $n $p $sigma $k $eta $alg $reps" run_simu.R $fnm &
done
done
done
done
done
done
done
