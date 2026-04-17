#!/bin/bash
set -eu

echo
echo ">>> TEST SMIFFER 6: Probe-Probe (PP) Fields"

fout="testdata/smiffer/pp_fields"
fpdb_nosolv="testdata/smiffer/pdb_clean"

mkdir -p $fout

conf_pp_base="DO_SMIF_APBS=False GRID_FORMAT_OUTPUT=MRC"

echo ">>> Testing protein PP fields generation..."

python3 volgrids smiffer prot $fpdb_nosolv/1iqj.pdb -o $fout -pp -c "$conf_pp_base"
echo "... PP fields generated"

# Verify PP field files were created
echo ">>> Verifying protein PP field outputs..."

if [ -f "$fout/1iqj.hbaPP.mrc" ]; then
    echo "... hbaPP field: OK"
else
    echo "... hbaPP field: MISSING"
fi

if [ -f "$fout/1iqj.hbdPP.mrc" ]; then
    echo "... hbdPP field: OK"
else
    echo "... hbdPP field: MISSING"
fi

if [ -f "$fout/1iqj.stkPP.mrc" ]; then
    echo "... stkPP field: OK"
else
    echo "... stkPP field: MISSING"
fi

if [ -f "$fout/1iqj.hpPP.mrc" ]; then
    echo "... hpPP field: OK"
else
    echo "... hpPP field: MISSING"
fi

############################# RNA PP FIELDS
echo ">>> Testing RNA PP fields generation..."

# Test basic PP fields for RNA
python3 volgrids smiffer rna $fpdb_nosolv/2esj.pdb -o $fout -pp -c "$conf_pp_base"
echo "... RNA PP fields generated"

# Verify RNA PP field files were created
echo ">>> Verifying RNA PP field outputs..."

if [ -f "$fout/2esj.hbaPP.mrc" ]; then
    echo "... RNA hbaPP field: OK"
else
    echo "... RNA hbaPP field: MISSING"
fi

if [ -f "$fout/2esj.hbdPP.mrc" ]; then
    echo "... RNA hbdPP field: OK"
else
    echo "... RNA hbdPP field: MISSING"
fi

if [ -f "$fout/2esj.stkPP.mrc" ]; then
    echo "... RNA stkPP field: OK"
else
    echo "... RNA stkPP field: MISSING"
fi

if [ -f "$fout/2esj.hpPP.mrc" ]; then
    echo "... RNA hpPP field: OK"
else
    echo "... RNA hpPP field: MISSING"
fi

############################# POCKET SPHERE PP FIELDS
echo ">>> Testing PP fields with pocket sphere..."

# Test PP fields with sphere constraint for protein
python3 volgrids smiffer prot $fpdb_nosolv/1iqj.pdb -o $fout -pp -s 4.682 21.475 7.161 14.675 -c "$conf_pp_base"
echo "... pocket sphere PP fields generated"

############################# PP FIELDS SUMMARIES
echo ">>> Grid summaries for PP fields..."

# Summarize protein PP fields
if [ -f "$fout/1iqj.hbaPP.mrc" ]; then
    python3 volgrids vgtools summary "$fout/1iqj.hbaPP.mrc"
fi

if [ -f "$fout/1iqj.hbdPP.mrc" ]; then
    python3 volgrids vgtools summary "$fout/1iqj.hbdPP.mrc"
fi

if [ -f "$fout/1iqj.stkPP.mrc" ]; then
    python3 volgrids vgtools summary "$fout/1iqj.stkPP.mrc"
fi

if [ -f "$fout/1iqj.hpPP.mrc" ]; then
    python3 volgrids vgtools summary "$fout/1iqj.hpPP.mrc"
fi

############################# COMPARISON WITH REGULAR FIELDS
echo ">>> Testing PP fields alongside regular fields..."

# Generate both regular and PP fields for comparison
python3 volgrids smiffer prot $fpdb_nosolv/1iqj.pdb -o $fout -pp -c "DO_SMIF_APBS=False DO_SMIF_HBA=True DO_SMIF_HBD=True DO_SMIF_STACKING=True DO_SMIF_HYDROPHOBIC=True GRID_FORMAT_OUTPUT=MRC"
echo "... regular and PP fields generated together"

# Verify both types exist
echo ">>> Verifying coexistence of regular and PP fields..."

if [ -f "$fout/1iqj.hbacceptors.mrc" ] && [ -f "$fout/1iqj.hbaPP.mrc" ]; then
    echo "... regular and PP HBA fields: OK"
else
    echo "... regular and PP HBA fields: ISSUE"
fi

if [ -f "$fout/1iqj.hbdonors.mrc" ] && [ -f "$fout/1iqj.hbdPP.mrc" ]; then
    echo "... regular and PP HBD fields: OK"
else
    echo "... regular and PP HBD fields: ISSUE"
fi
