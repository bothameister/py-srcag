#!/bin/bash

for f in `ls *.yld`; do
  base=${f%%.yld}
  rm -f $base.g
  ../preprocess-grammar.py < $base > $base.g
  echo -e "\n===== $base ====="
  ../py-srcag -d 10000 -q $base.g < $f
done

rm *.g

