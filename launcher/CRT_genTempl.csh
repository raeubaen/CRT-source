#!/bin/tcsh

cd ../..

g++ ./src/CRT_genTemplate.C `root-config --cflags --glibs` -o ./x/CRT_genTemplate.x

./x/CRT_genTemplate.x 0 0 0 0

