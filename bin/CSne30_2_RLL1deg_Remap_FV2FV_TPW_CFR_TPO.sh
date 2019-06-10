# Generate the mapping data for testing CSne30 to RLL1deg
source envs.sh
$TEMPESTREMAP_DIR/./GenerateOverlapMesh --a ../meshes/outCSne30.g --b ../meshes/outRLL1deg.g --out ../meshes/CSne30_2_RLL1deg.g
$TEMPESTREMAP_DIR/./GenerateOfflineMap --in_mesh ../meshes/outCSne30.g --out_mesh ../meshes/outRLL1deg.g --ov_mesh ../meshes/CSne30_2_RLL1deg.g --in_np 4 --out_map ../meshes/CSne30_2_RLL1deg_np4.nc

# Remap TPW from CSne30 to RLL1deg
$TEMPESTREMAP_DIR/./ApplyOfflineMap --map ../meshes/CSne30_2_RLL1deg_np4.nc --var TotalPrecipWater --in_data ../testdata_outCSne30_TPW_CFR_TPO.nc --out_data ../testdata_outCSne30_2_RLL1deg_TPW.nc
# Remap CFR from CSne30 to RLL1deg
$TEMPESTREMAP_DIR/./ApplyOfflineMap --map ../meshes/CSne30_2_RLL1deg_np4.nc --var CloudFraction --in_data ../testdata_outCSne30_TPW_CFR_TPO.nc --out_data ../testdata_outCSne30_2_RLL1deg_CFR.nc
# Remap TPO from CSne30 to RLL1deg
$TEMPESTREMAP_DIR/./ApplyOfflineMap --map ../meshes/CSne30_2_RLL1deg_np4.nc --var Topography --in_data ../testdata_outCSne30_TPW_CFR_TPO.nc --out_data ../testdata_outCSne30_2_RLL1deg_TPO.nc
