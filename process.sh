#!/bin/bash
pdb_code=$1
if test ! -e $pdb_code.pdb; then
phenix.fetch_pdb $pdb_code > /dev/null
fi
if test ! -e $pdb_code-sf.cif; then
phenix.fetch_pdb -x $pdb_code > /dev/null
fi
if test ! -e $pdb_code.ligands.cif; then
phenix.ready_set optimise_final_geometry_of_hydrogens=False $pdb_code.pdb > /dev/null
fi
threshold=4.0
rmsd=`phenix.pdbtools model_statistics=true stop_for_unknowns=false $pdb_code.pdb | grep Mean | awk '{print $2}'`
echo $pdb_code' '$rmsd >> ../rmsd_log
if (( `echo $rmsd'>'$threshold | bc -l` )); then
	echo $pdb_code >> ../b_sharp
	phenix.maps $pdb_code-sf.cif $pdb_code.pdb ../maps.params | grep b_sharp: >> ../b_sharp
	if test -e $pdb_code.ligands.cif; then
		phenix.get_cc_mtz_pdb $pdb_code\_map_coeffs.mtz $pdb_code.pdb $pdb_code.ligands.cif \
			labin="FP=FWT PHIB=PHFWT" fix_xyz=True | grep CC
		phenix.get_cc_mtz_pdb $pdb_code\_map_coeffs.mtz $pdb_code.pdb $pdb_code.ligands.cif \
			labin="FP=2FOFCWT PHIB=PH2FOFCWT" fix_xyz=True | grep CC
	else
		phenix.get_cc_mtz_pdb $pdb_code\_map_coeffs.mtz $pdb_code.pdb \
			labin="FP=FWT PHIB=PHFWT" fix_xyz=True | grep CC
		phenix.get_cc_mtz_pdb $pdb_code\_map_coeffs.mtz $pdb_code.pdb  \
			labin="FP=2FOFCWT PHIB=PH2FOFCWT" fix_xyz=True | grep CC
	fi
fi
