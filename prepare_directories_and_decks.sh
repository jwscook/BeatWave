#!/bin/bash

datestub=$(date +%Y-%m-%d-%H-%M-%S)

mkdir -p $datestub
mkdir -p $datestub/null
mkdir -p $datestub/wave1
mkdir -p $datestub/wave2
mkdir -p $datestub/both
for i in null wave1 wave2 both
do
  cat $i.prefix.deck suffix.deck > $datestub/$i/input.deck
done

# Run via e.g.
# echo 2023-01-04-09-39-23/both/ | mpirun -np 2 ../epoch/epoch1d/bin/epoch1d
