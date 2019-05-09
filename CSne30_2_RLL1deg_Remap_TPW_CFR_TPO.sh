# Generate the mapping data for testing CSne30 to RLL1deg
~/Desktop/tempestremap/build/./GenerateOverlapMesh --a outCSne30.g --b outRLL1deg.g --out CSne30_2_RLL1deg.g
~/Desktop/tempestremap/build/./GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outRLL1deg.g --ov_mesh CSne30_2_RLL1deg.g --in_np 4 --out_map CSne30_2_RLL1deg_np4_TPW_CFR_TPO.nc

# Remap TPW from CSne30 to RLL1deg
~/Desktop/tempestremap/build/./ApplyOfflineMap --map CSne30_2_RLL1deg_np4_TPW_CFR_TPO.nc --var TotalPrecipWater --in_data testdata_outCSne30_TPW_CFR_TPO.nc --out_data testdata_outCSne30_2_RLL1deg_TPW.nc
# Remap CFR from CSne30 to RLL1deg
~/Desktop/tempestremap/build/./ApplyOfflineMap --map CSne30_2_RLL1deg_np4_TPW_CFR_TPO.nc --var CloudFraction --in_data testdata_outCSne30_TPW_CFR_TPO.nc --out_data testdata_outCSne30_2_RLL1deg_CFR.nc
# Remap TPO from CSne30 to RLL1deg
~/Desktop/tempestremap/build/./ApplyOfflineMap --map CSne30_2_RLL1deg_np4_TPW_CFR_TPO.nc --var Topography --in_data testdata_outCSne30_TPW_CFR_TPO.nc --out_data testdata_outCSne30_2_RLL1deg_TPO.nc
