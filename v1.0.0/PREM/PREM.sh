#!/bin/bash

obj="PREM"
rm -f $obj
g++ -std=c++11 -O3 $obj.cxx -o $obj
time ./$obj
rm -f $obj
