# Generate the mapping data for testing CSne60 to ICO64
source envs.sh
$TEMPESTREMAP_DIR/GenerateOverlapMesh --a outCSne60.g --b outICO64.g --out CSne60_2_ICO64.g
$TEMPESTREMAP_DIR/GenerateOfflineMap --in_mesh outCSne60.g --out_mesh outICO64.g --ov_mesh CSne60_2_ICO64.g --in_np 4 --out_map CSne60_2_ICO64_np4_TPW_CFR_TPO.nc

# Remap TPW from CSne60 to ICO64
$TEMPESTREMAP_DIR/ApplyOfflineMap --map CSne60_2_ICO64_np4_TPW_CFR_TPO.nc --var TotalPrecipWater --in_data testdata_outCSne60_TPW_CFR_TPO.nc --out_data testdata_outCSne60_2_ICO64_TPW.nc
# Remap CFR from CSne60 to ICO64
$TEMPESTREMAP_DIR/ApplyOfflineMap --map CSne60_2_ICO64_np4_TPW_CFR_TPO.nc --var CloudFraction --in_data testdata_outCSne60_TPW_CFR_TPO.nc --out_data testdata_outCSne60_2_ICO64_CFR.nc
# Remap TPO from CSne60 to ICO64
$TEMPESTREMAP_DIR/ApplyOfflineMap --map CSne60_2_ICO64_np4_TPW_CFR_TPO.nc --var Topography --in_data testdata_outCSne60_TPW_CFR_TPO.nc --out_data testdata_outCSne60_2_ICO64_TPO.nc
