#!/usr/bin/env bash

START=$(date +%s.%N)
python fftdemo.py
END=$(date +%s.%N)
DIFF=$(echo "$END - $START" | bc)

echo "dft time = ${DIFF}"

START=$(date +%s.%N)
python fftdemo2.py
END=$(date +%s.%N)
DIFF=$(echo "$END - $START" | bc)

echo "fft time = $DIFF"
