#!/bin/bash
while read -r line || [[ -n "$line" ]]; do
  if [ -d "$line" ]; then
  	cd $line
	{ echo -n "$line " ; phenix.pdbtools model_statistics=true stop_for_unknowns=false $line.pdb | grep Mean ; } >> ../stat_log
  	cd ..
  fi
done < "$1"
