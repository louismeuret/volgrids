#!/bin/bash
set -eu

##### Generate input data (if needed)
if [[ ! -f testdata/smiffer/apbs/1akx.pdb.mrc ]]; then
    echo "Generating input data..."
    bash tests/_gen_input.sh
else
    echo "Input data already exists, skipping generation."
fi


##### Run tests
bash tests/env/vgtest.sh

bash tests/smiffer/toy_systems.sh
bash tests/smiffer/pocket_sphere.sh
bash tests/smiffer/whole.sh
bash tests/smiffer/traj.sh
bash tests/smiffer/ligand.sh
bash tests/smiffer/cavities.sh
bash tests/smiffer/pp_fields.sh

bash tests/vgtools/convert.sh
bash tests/vgtools/pack_unpack.sh
bash tests/vgtools/fix_cmap.sh
bash tests/vgtools/compare.sh
bash tests/vgtools/overlap.sh

echo "All tests completed successfully."
