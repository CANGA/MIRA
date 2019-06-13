# Generate the mapping data for testing CSne60 to RLL1deg in FV2FV
source envs.sh
$TEMPESTREMAP_DIR/./GenerateOverlapMesh --a ../meshes/outCSne60.g --b ../meshes/outRLL1deg.g --out ../meshes/CSne60_2_RLL1deg.g
$TEMPESTREMAP_DIR/./GenerateOfflineMap --in_mesh ../meshes/outCSne60.g --out_mesh ../meshes/outRLL1deg.g --ov_mesh ../meshes/CSne60_2_RLL1deg.g --in_np 4 --out_map ../meshes/CSne60_2_RLL1deg.nc

# Remap test 3 from CSne60 to RLL1deg
$TEMPESTREMAP_DIR/./ApplyOfflineMap --map ../meshes/CSne60_2_RLL1deg.nc --var Psi --in_data ../testdata_CSne60_np4_3.nc --out_data ../testdata_outCSne60_2_RLL1deg_np4_3.nc
