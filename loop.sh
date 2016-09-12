#!/bin/bash
while read -r line || [[ -n "$line" ]]; do
  mkdir -p $line
  cd $line
  ../process.sh $line 1>> ../log 2>> ../err
  cd ..
done < "$1"
