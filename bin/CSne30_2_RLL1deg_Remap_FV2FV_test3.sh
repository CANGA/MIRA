# Generate the mapping data for testing CSne30 to RLL1deg in FV2FV
source envs.sh
$TEMPESTREMAP_DIR/./GenerateOverlapMesh --a ../meshes/outCSne30.g --b ../meshes/outRLL1deg.g --out ../meshes/CSne30_2_RLL1deg.g
$TEMPESTREMAP_DIR/./GenerateOfflineMap --in_mesh ../meshes/outCSne30.g --out_mesh ../meshes/outRLL1deg.g --ov_mesh ../meshes/CSne30_2_RLL1deg.g --in_np 4 --out_map ../meshes/CSne30_2_RLL1deg.nc

# Remap test 3 from CSne30 to RLL1deg
$TEMPESTREMAP_DIR/./ApplyOfflineMap --map ../meshes/CSne30_2_RLL1deg.nc --var Psi --in_data ../testdata_CSne30_np4_3.nc --out_data ../testdata_outCSne30_2_RLL1deg_np4_3.nc
