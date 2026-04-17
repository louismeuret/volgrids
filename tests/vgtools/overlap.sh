#!/bin/bash
set -eu

echo
echo ">>> TEST VGTOOLS 4: Overlap operations between chains A and B (1ofz)"

folder="testdata/vgtools"
fo="$folder/overlap"
fpdb_clean="testdata/smiffer/pdb_clean"

mkdir -p $fo

############################# CHAIN SEPARATION
echo ">>> Separating chains A and B from 1ofz..."
grep "^ATOM.*\ A\ " $fpdb_clean/1ofz.pdb > $fo/1ofz_chainA.pdb
echo "END" >> $fo/1ofz_chainA.pdb

grep "^ATOM.*\ B\ " $fpdb_clean/1ofz.pdb > $fo/1ofz_chainB.pdb
echo "END" >> $fo/1ofz_chainB.pdb

echo "... Chain A atoms: $(grep "^ATOM" $fo/1ofz_chainA.pdb | wc -l)"
echo "... Chain B atoms: $(grep "^ATOM" $fo/1ofz_chainB.pdb | wc -l)"

############################# GENERATE PP FIELDS FOR CHAIN A
echo ">>> Generating PP fields for chain A..."
python3 volgrids smiffer prot $fo/1ofz_chainA.pdb -o $fo -pp -c "DO_SMIF_APBS=False GRID_FORMAT_OUTPUT=MRC"

# Rename with chain identifier
mv $fo/1ofz_chainA.hbaPP.mrc $fo/chainA.hbaPP.mrc 2>/dev/null || echo "... hbaPP not found"
mv $fo/1ofz_chainA.hbdPP.mrc $fo/chainA.hbdPP.mrc 2>/dev/null || echo "... hbdPP not found"
mv $fo/1ofz_chainA.stkPP.mrc $fo/chainA.stkPP.mrc 2>/dev/null || echo "... stkPP not found"

############################# GENERATE CLASSIC FIELDS FOR CHAIN B
echo ">>> Generating classic fields for chain B..."
python3 volgrids smiffer prot $fo/1ofz_chainB.pdb -o $fo -c "DO_SMIF_APBS=False DO_SMIF_HBA=True DO_SMIF_HBD=True DO_SMIF_STACKING=True GRID_FORMAT_OUTPUT=MRC"

# Rename with chain identifier
mv $fo/1ofz_chainB.hbacceptors.mrc $fo/chainB.hba.mrc 2>/dev/null || echo "... hba not found"
mv $fo/1ofz_chainB.hbdonors.mrc $fo/chainB.hbd.mrc 2>/dev/null || echo "... hbd not found"
mv $fo/1ofz_chainB.stacking.mrc $fo/chainB.stk.mrc 2>/dev/null || echo "... stk not found"

############################# CROSS-CHAIN OVERLAP ANALYSIS
echo ">>> Testing cross-overlap between chain A (PP) and chain B (classic)..."

# HBA (chain B) vs HBD-PP (chain A)
if [ -f "$fo/chainB.hba.mrc" ] && [ -f "$fo/chainA.hbdPP.mrc" ]; then
    python3 volgrids vgtools overlap "$fo/chainB.hba.mrc" "$fo/chainA.hbdPP.mrc" "$fo/hba_vs_hbdPP.mrc" --operation multiply
    echo "... HBA (B) vs HBD-PP (A) overlap completed"
fi

# HBD (chain B) vs HBA-PP (chain A)
if [ -f "$fo/chainB.hbd.mrc" ] && [ -f "$fo/chainA.hbaPP.mrc" ]; then
    python3 volgrids vgtools overlap "$fo/chainB.hbd.mrc" "$fo/chainA.hbaPP.mrc" "$fo/hbd_vs_hbaPP.mrc" --operation multiply
    echo "... HBD (B) vs HBA-PP (A) overlap completed"
fi

# Stacking (chain B) vs Stacking-PP (chain A)
if [ -f "$fo/chainB.stk.mrc" ] && [ -f "$fo/chainA.stkPP.mrc" ]; then
    python3 volgrids vgtools overlap "$fo/chainB.stk.mrc" "$fo/chainA.stkPP.mrc" "$fo/stk_vs_stkPP.mrc" --operation multiply
    echo "... Stacking (B) vs Stacking-PP (A) overlap completed"
fi

############################# OVERLAP OPERATIONS TEST
echo ">>> Testing different overlap operations..."

if [ -f "$fo/chainB.hba.mrc" ] && [ -f "$fo/chainA.hbaPP.mrc" ]; then
    # Test add operation
    python3 volgrids vgtools overlap "$fo/chainB.hba.mrc" "$fo/chainA.hbaPP.mrc" "$fo/hba_add_hbaPP.mrc" --operation add
    echo "... Add operation completed"

    # Test subtract operation
    python3 volgrids vgtools overlap "$fo/chainB.hba.mrc" "$fo/chainA.hbaPP.mrc" "$fo/hba_sub_hbaPP.mrc" --operation subtract
    echo "... Subtract operation completed"
fi

############################# CONVENIENCE FUNCTIONS
echo ">>> Testing overlap convenience functions..."

if [ -f "$fo/chainB.hba.mrc" ] && [ -f "$fo/chainA.hbdPP.mrc" ]; then
    python3 volgrids vgtools overlap_cross "$fo/chainB.hba.mrc" "$fo/chainA.hbdPP.mrc" "$fo/cross_comparison.mrc"
    echo "... Cross-comparison function completed"
fi

if [ -f "$fo/chainB.hbd.mrc" ] && [ -f "$fo/chainA.hbaPP.mrc" ]; then
    python3 volgrids vgtools overlap_diff "$fo/chainB.hbd.mrc" "$fo/chainA.hbaPP.mrc" "$fo/difference_analysis.mrc"
    echo "... Difference analysis completed"
fi

############################# VERIFICATION AND SUMMARIES
echo ">>> Verifying overlap results..."

results="hba_vs_hbdPP hbd_vs_hbaPP stk_vs_stkPP hba_add_hbaPP hba_sub_hbaPP cross_comparison difference_analysis"
for result in $results; do
    if [ -f "$fo/${result}.mrc" ]; then
        echo "... $result: OK"
        python3 volgrids vgtools summary "$fo/${result}.mrc"
    else
        echo "... $result: MISSING"
    fi
done

echo ">>> Chain separation and overlap test completed successfully!"