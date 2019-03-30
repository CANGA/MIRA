# Run metrics from CSne30 to RLL1deg on Total Precip. Water, Cloud Fraction, and Terrain
python CANGAMetricsDriver.py -v TotalPrecipWater --ss testdata_outCSne30_TPW_CFR_TPO.nc --s2t testdata_outCSne30_2_RLL1deg_TPW.nc --st testdata_outRLL1deg_TPW_CFR_TPO.nc --sm outCSne30.g --tm outRLL1deg.g --ExodusSingleConn

python CANGAMetricsDriver.py -v CloudFraction --ss testdata_outCSne30_TPW_CFR_TPO.nc --s2t testdata_outCSne30_2_RLL1deg_CFR.nc --st testdata_outRLL1deg_TPW_CFR_TPO.nc --sm outCSne30.g --tm outRLL1deg.g --ExodusSingleConn

python CANGAMetricsDriver.py -v Topography --ss testdata_outCSne30_TPW_CFR_TPO.nc --s2t testdata_outCSne30_2_RLL1deg_TPO.nc --st testdata_outRLL1deg_TPW_CFR_TPO.nc --sm outCSne30.g --tm outRLL1deg.g --ExodusSingleConn
