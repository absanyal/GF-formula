#!/bin/bash
rm -r w*
for i in {1..100}
do
    WMIN=-5
    WMAX=45
    DELTA1=$((WMAX-WMIN))
    DELTA=$(div $DELTA1 100)
    W=`echo $WMIN+$DELTA\*$i|bc`
    # echo $DELTA
    echo $W
    mkdir w${i}
    cd w${i}
    echo $W > omega_value.dat
    cp ../spinless.py spinless.py
    cp ../statemanip.py statemanip.py
    cp ../job.sh job.sh
    chmod 777 job.sh
    qsub job.sh
    cd ../
done