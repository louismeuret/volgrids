#!/bin/bash
set -eu

echo
echo ">>> TEST SMUTILS 0: Watson-Crick base pair detection"

python3 - <<'EOF'
from volgrids.smutils._core.rna_resids import RNAResids
from pathlib import Path

pdb = Path("testdata/smiffer/whole/1akx.pdb")
result = RNAResids.get_resids_bp_wc(pdb)

assert isinstance(result, set), "Expected a set"
assert len(result) > 0, "Expected at least one WC residue"
assert all(isinstance(r, int) for r in result), "Expected integer resids"

print(f"  WC resids ({len(result)}): {sorted(result)}")
print("  PASS")
EOF
