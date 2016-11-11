#!/bin/bash

gcc $3 -lpthread -lrt -lm

for(( i = $1; i < $2; i++))
do

	./a.out $i

done
