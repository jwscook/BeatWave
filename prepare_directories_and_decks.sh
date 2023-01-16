#!/bin/bash

if [ -z "$1" ]
  then
    dirname=$(date +%Y-%m-%d-%H-%M-%S)
else
    dirname=$1
fi

mkdir -p $dirname
mkdir -p $dirname/null
mkdir -p $dirname/wave1
mkdir -p $dirname/wave2
mkdir -p $dirname/both
for i in null wave1 wave2 both
do
  cat $i.prefix.deck suffix.deck > $dirname/$i/input.deck
done

# Run via e.g.
# echo 2023-01-04-09-39-23/both/ | mpirun -np 2 ../epoch/epoch1d/bin/epoch1d
