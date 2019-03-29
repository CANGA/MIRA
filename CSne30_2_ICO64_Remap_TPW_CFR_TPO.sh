# Generate the mapping data for testing CSne30 to ICO64
~/tempestremap/bin/./GenerateOverlapMesh --a outCSne30.g --b outICO64.g --out CSne30_2_ICO64.g
~/tempestremap/bin/./GenerateOfflineMap --in_mesh outCSne30.g --out_mesh outICO64.g --ov_mesh CSne30_2_ICO64.g --in_np 4 --out_map CSne30_2_ICO64_np4_TPW_CFR_TPO.nc

# Remap TPW from CSne30 to ICO64
~/tempestremap/bin/./ApplyOfflineMap --map CSne30_2_ICO64_np4_TPW_CFR_TPO.nc --var TotalPrecipWater --in_data testdata_outCSne30_TPW_CFR_TPO.nc --out_data testdata_outCSne30_2_ICO64_TPW.nc
# Remap TPW from CSne30 to ICO64
~/tempestremap/bin/./ApplyOfflineMap --map CSne30_2_ICO64_np4_TPW_CFR_TPO.nc --var CloudFraction --in_data testdata_outCSne30_TPW_CFR_TPO.nc --out_data testdata_outCSne30_2_ICO64_CFR.nc
# Remap TPW from CSne30 to ICO64
~/tempestremap/bin/./ApplyOfflineMap --map CSne30_2_ICO64_np4_TPW_CFR_TPO.nc --var Topography --in_data testdata_outCSne30_TPW_CFR_TPO.nc --out_data testdata_outCSne30_2_ICO64_TPO.nc
