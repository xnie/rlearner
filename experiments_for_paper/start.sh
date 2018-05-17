#! /bin/bash

setupvals=(1 2 3 4 5 6 7 8)
nvals=(500 1000)
pvals=(6 12)
sigmavals=(0.5 1.0 2.0 4.0)
algvals=('R' 'X' 'S' 'T' 'U' 'RS' 'oracle')
lambdachoice='lambda.min'

reps=200

for ((i1=0; i1<${#setupvals[@]} ;i1++))
do
for ((i2=0; i2<${#nvals[@]} ;i2++))
do
for ((i3=0; i3<${#pvals[@]} ;i3++))
do
for ((i4=0; i4<${#sigmavals[@]} ;i4++))
do
for ((i5=0; i5<${#algvals[@]} ;i5++))
do
    while [ `pgrep -c R` -ge 100 ]
    do
          sleep 5
    done

    setup=${setupvals[$i1]}
    n=${nvals[$i2]}
    p=${pvals[$i3]}
    sigma=${sigmavals[$i4]}
    alg=${algvals[$i5]}

    fnm="logging/progress-$alg-$setup-$n-$p-$sigma-$reps-$lambdachoice.out"
    #echo $fnm

    #Rscript run_simu.R $alg $setup $n $p $sigma $reps $lambdachoice 2>&1 | tee $fnm &
    echo "Rscript run_simu.R $alg $setup $n $p $sigma $reps $lambdachoice 2>&1 | tee $fnm &"
    #R CMD BATCH --no-save --no-restore "--args $setup $n $p $sigma $k $eta $alg $reps" run_simu.R $fnm &
done
done
done
done
done
