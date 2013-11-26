#!/usr/bin/env bash

i=10
while [ $i -lt 10000 ]
do
    ./timeit.sh "$(echo "$i * 2" | bc)" "$(echo "$i * 10" | bc)"
    i="$(echo "$i + 10" | bc)"
done
