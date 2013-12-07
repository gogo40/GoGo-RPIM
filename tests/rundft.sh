#!/bin/bash
g++ -O2 -o dft dft.cpp

read n
read freq

for (( i = 0 ; i < 4 ; i++ ))
do
	read fname
	./dft $freq $fname
	mv F.dat.dft F.$i.dat
done

