# Generate the mapping data for testing RLL1deg to CSne30 SE mode
source envs.sh
$TEMPESTREMAP_DIR/./GenerateOverlapMesh --b ../meshes/outCSne30.g --a ../meshes/outRLL1deg.g --out ../meshes/RLL1deg_2_CSne30.g
$TEMPESTREMAP_DIR/./GenerateOfflineMap --out_mesh ../meshes/outCSne30.g --in_mesh ../meshes/outRLL1deg.g --ov_mesh ../meshes/RLL1deg_2_CSne30.g --in_np 4 --in_type fv --out_type cgll --out_map ../meshes/RLL1deg_2_CSne30_np4.nc

# Remap TPW from CSne30 to RLL1deg
$TEMPESTREMAP_DIR/./ApplyOfflineMap --map ../meshes/RLL1deg_2_CSne30_np4.nc --var TotalPrecipWater --in_data ../testdata_outRLL1deg_TPW_CFR_TPO.nc --out_data ../testdataGLL_outCSne30_2_RLL1deg_TPW.nc
# Remap CFR from CSne30 to RLL1deg
$TEMPESTREMAP_DIR/./ApplyOfflineMap --map ../meshes/RLL1deg_2_CSne30_np4.nc --var CloudFraction --in_data ../testdata_outRLL1deg_TPW_CFR_TPO.nc --out_data ../testdataGLL_outCSne30_2_RLL1deg_CFR.nc
# Remap TPO from CSne30 to RLL1deg
$TEMPESTREMAP_DIR/./ApplyOfflineMap --map ../meshes/RLL1deg_2_CSne30_np4.nc --var Topography --in_data ../testdata_outRLL1deg_TPW_CFR_TPO.nc --out_data ../testdataGLL_outCSne30_2_RLL1deg_TPO.nc
