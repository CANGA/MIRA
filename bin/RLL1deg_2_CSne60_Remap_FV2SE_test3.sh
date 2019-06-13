# Generate the mapping data for testing RLL1deg to CSne60 SE mode
source envs.sh
$TEMPESTREMAP_DIR./GenerateOverlapMesh --b ../meshes/outCSne60.g --a ../meshes/outRLL1deg.g --out ../meshes/RLL1deg_2_CSne60.g
$TEMPESTREMAP_DIR./GenerateOfflineMap --in_mesh ../meshes/outRLL1deg.g --out_mesh ../meshes/outCSne60.g --ov_mesh ../meshes/RLL1deg_2_CSne60.g --in_type fv --out_type cgll --in_np 4 --out_np 4 --out_map ../meshes/RLL1deg_2_CSne60_np4.nc

# Remap test 3 from RLL1deg to CSne60
$TEMPESTREMAP_DIR./ApplyOfflineMap --map ../meshes/RLL1deg_2_CSne60_np4.nc --var Psi --in_data ../testdata_RLL1deg_np4_3.nc --out_data ../testdataGLL_RLL1deg_2_outCSne60_np4_3.nc
