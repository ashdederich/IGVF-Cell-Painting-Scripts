#!/bin/bash
module add python/3.8.x-anaconda
conda activate useful_packages

dir=$1

cd $dir

find . -name "* *" -type f | rename 's/ /_/g'
for f in *; do mv "$f" "${f//\(/_}"; done
for f in *; do mv "$f" "${f//\)/}"; done

cd ..
