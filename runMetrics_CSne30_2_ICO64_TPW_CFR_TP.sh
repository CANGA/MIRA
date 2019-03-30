# Run metrics from CSne30 to ICO64 on Total Precip. Water, Cloud Fraction, and Terrain
python CANGAMetricsDriver.py -v TotalPrecipWater --ss testdata_outCSne30_TPW_CFR_TPO.nc --s2t testdata_outCSne30_2_ICO64_TPW.nc --st testdata_outICO64_TPW_CFR_TPO.nc --sm outCSne30.g --tm outICO64.g --ExodusSingleConn

python CANGAMetricsDriver.py -v CloudFraction --ss testdata_outCSne30_TPW_CFR_TPO.nc --s2t testdata_outCSne30_2_ICO64_CFR.nc --st testdata_outICO64_TPW_CFR_TPO.nc --sm outCSne30.g --tm outICO64.g --ExodusSingleConn

python CANGAMetricsDriver.py -v Topography --ss testdata_outCSne30_TPW_CFR_TPO.nc --s2t testdata_outCSne30_2_ICO64_TPO.nc --st testdata_outICO64_TPW_CFR_TPO.nc --sm outCSne30.g --tm outICO64.g --ExodusSingleConn
