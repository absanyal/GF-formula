#!/bin/bash
rm -r w*
for i in {1..100}
do
    WMIN=-5
    WMAX=45
    DELTA1=$((WMAX-WMIN))
    DELTA=`echo $DELTA1/100|bc -l`
    W=`echo $WMIN+$DELTA\*$i|bc`
    # echo $DELTA
    # echo $W
    mkdir w${i}
    cd w${i}
    echo $W > omega_value.txt
    cp ../spinless.py spinless.py
    cp ../statemanip.py statemanip.py
    cp ../submit.sh submit.sh
    chmod 777 submit.sh
    qsub submit.sh
    # python3 spinless.py
    cd ../
done