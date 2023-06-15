#!/bin/bash

n1=`expr ${1} - 1`
n2=`expr ${2} - 1`
n3=`expr ${3} - 1`

echo "K_POINTS crystal"
expr ${1} \* ${2} \* ${3} \* 2

for i1 in `seq 0 ${n1}`; do
    kv1=`echo "scale = 15; ${i1} / ${1}" | bc`
    if [ `expr ${i1} \* 2` -ge ${1} ]; then
        kv1=`echo "scale = 15; ${kv1} - 1.0" | bc`
    fi
    for i2 in `seq 0 ${n2}`; do
        kv2=`echo "scale = 15; ${i2} / ${2}" | bc`
        if [ `expr ${i2} \* 2` -ge ${2} ]; then
            kv2=`echo "scale = 15; ${kv2} - 1.0" | bc`
        fi
        for i3 in `seq 0 ${n3}`; do
            kv3=`echo "scale = 15; ${i3} / ${3}" | bc`
            if [ `expr ${i3} \* 2` -ge ${3} ]; then
                kv3=`echo "scale = 15; ${kv3} - 1.0" | bc`
            fi

            printf "%20.15f %20.15f %20.15f   1.0\n" ${kv1} ${kv2} ${kv3}
        done
    done
done

for i1 in `seq 0 ${n1}`; do
    kv1=`echo "scale = 15; (2 * ${i1} + 1) / (2 * ${1})" | bc`
    if [ `expr ${i1} \* 2 + 1` -ge ${1} ]; then
        kv1=`echo "scale = 15; ${kv1} - 1.0" | bc`
    fi
    for i2 in `seq 0 ${n2}`; do
        kv2=`echo "scale = 15; (2 * ${i2} + 1) / (2 * ${2})" | bc`
        if [ `expr ${i2} \* 2 + 1` -ge ${2} ]; then
            kv2=`echo "scale = 15; ${kv2} - 1.0" | bc`
        fi
        for i3 in `seq 0 ${n3}`; do
            kv3=`echo "scale = 15; (2 * ${i3} + 1) / (2 * ${3})" | bc`
            if [ `expr ${i3} \* 2 + 1` -ge ${3} ]; then
                kv3=`echo "scale = 15; ${kv3} - 1.0" | bc`
            fi
            
            printf "%20.15f %20.15f %20.15f   1.0\n" ${kv1} ${kv2} ${kv3}
        done
    done
done
