#! /bin/bash

setupvals=('A' 'B' 'C' 'D' 'E' 'F')
nvals=(500 1000)
pvals=(6 12)
sigmavals=(0.5 1.0 2.0 4.0)
lassoalgvals=('R' 'RS' 'T' 'X' 'U' 'oracle' 'S')
boostalgvals=('R' 'RC' 'T' 'X' 'U' 'oracle' 'S' 'causalboost')

reps=100
learners=('boost' 'lasso')

for ((i1=0; i1<${#setupvals[@]} ;i1++))
do
for ((i2=0; i2<${#nvals[@]} ;i2++))
do
for ((i3=0; i3<${#pvals[@]} ;i3++))
do
for ((i4=0; i4<${#sigmavals[@]} ;i4++))
do
for ((i5=0; i5<${#learners[@]} ;i5++))
do
    while [ `pgrep -c R` -ge 100 ]
    do
        sleep 5
    done

    setup=${setupvals[$i1]}
    n=${nvals[$i2]}
    p=${pvals[$i3]}
    sigma=${sigmavals[$i4]}
    learner=${learners[$i5]}

    if [ "$learner" = "boost" ]; then
      algvals=( "${boostalgvals[@]}" )
    elif [ "$learner" = "lasso" ]; then
      algvals=( "${lassoalgvals[@]}" )
    else
      echo "learner needs to be lasso for boost for the experiments.";
      exit 1;
    fi

    for ((i6=0; i6<${#algvals[@]} ;i6++))
    do
      alg=${algvals[$i6]}
      fnm="logging/progress-$alg-$learner-$setup-$n-$p-$sigma-$reps.out"

      Rscript run_simu.R $alg $learner $setup $n $p $sigma $reps 2>&1 | tee $fnm &
      #echo "Rscript run_simu.R $alg $learner $setup $n $p $sigma $reps 2>&1 | tee $fnm &"
    done
done
done
done
done
done
