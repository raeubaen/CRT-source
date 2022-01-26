#!/bin/tcsh

cd ./../..

g++ ./src/CRT_step3p1.C `root-config --cflags --glibs` -o ./x/CRT_step3p1.x
g++ ./src/CRT_step3p2.C `root-config --cflags --glibs` -o ./x/CRT_step3p2.x

./x/CRT_step3p1.x ./data/step2/$1_s2.root ./data/step3p/$1_s3p1.root $1 $1
./x/CRT_step3p2.x ./data/step2/$1_s2.root ./data/step3p/$1_s3p2.root $1 $1

