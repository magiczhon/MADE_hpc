#!/bin/bash
make
list_N=( 500 512 1000 1024 2000 2048 )
echo > time.txt
for idx in `seq 0 $(( ${#list_N[@]}-1 ))`
do
	echo N = ${list_N[$idx]}
	echo ${list_N[$idx]} | ./result.exe >> time.txt
done