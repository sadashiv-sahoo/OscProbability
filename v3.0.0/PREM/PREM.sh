#!/bin/bash

Target="PREM"
rm -f $Target
g++ -std=c++11 -O3 $Target.cxx -o $Target
time ./$Target
rm -f $Target
