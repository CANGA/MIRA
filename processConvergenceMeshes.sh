# Run metrics from CSne30 to RLL1deg on Total Precip. Water, Cloud Fraction, and Terrain
python CANGAMeshPreprocessDriver.py --mesh convergenceMeshes/outCSne16.g --ExodusSingleConn
python CANGAMeshPreprocessDriver.py --mesh convergenceMeshes/outCSne32.g --ExodusSingleConn
python CANGAMeshPreprocessDriver.py --mesh convergenceMeshes/outCSne64.g --ExodusSingleConn
python CANGAMeshPreprocessDriver.py --mesh convergenceMeshes/outCSne128.g --ExodusSingleConn