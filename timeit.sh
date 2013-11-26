#!/usr/bin/env bash

echo "FFT Time elapsed: "
START=$(ruby -e 'puts "%.3f" % Time.now')
./fft $1 $2 > /dev/null
END=$(ruby -e 'puts "%.3f" % Time.now')
echo $(echo "$END - $START" | bc) seconds
echo "DFT Time elapsed: "
START=$(ruby -e 'puts "%.3f" % Time.now')
./dft $1 $2 > /dev/null
END=$(ruby -e 'puts "%.3f" % Time.now')
echo $(echo "$END - $START" | bc) seconds
