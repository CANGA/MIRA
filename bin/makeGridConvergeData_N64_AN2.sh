# Sample terrain with constant number of modes at integral order 4
python CANGAFieldGenDriver.py --pm convergenceMeshes/outCSne16.g --so 4 --nm 64 --evaluateA2 --ExodusSingleConn
python CANGAFieldGenDriver.py --pm convergenceMeshes/outCSne32.g --so 4 --nm 64 --evaluateA2 --ExodusSingleConn
python CANGAFieldGenDriver.py --pm convergenceMeshes/outCSne64.g --so 4 --nm 64 --evaluateA2 --ExodusSingleConn
python CANGAFieldGenDriver.py --pm convergenceMeshes/outCSne128.g --so 4 --nm 64 --evaluateA2 --ExodusSingleConn
# Sample terrain with constant number of modes at integral order 1 (centroid sampling)
python CANGAFieldGenDriver.py --pm convergenceMeshes/outCSne16.g --so 1 --nm 64 --evaluateA2 --ExodusSingleConn
python CANGAFieldGenDriver.py --pm convergenceMeshes/outCSne32.g --so 1 --nm 64 --evaluateA2 --ExodusSingleConn
python CANGAFieldGenDriver.py --pm convergenceMeshes/outCSne64.g --so 1 --nm 64 --evaluateA2 --ExodusSingleConn
python CANGAFieldGenDriver.py --pm convergenceMeshes/outCSne128.g --so 1 --nm 64 --evaluateA2 --ExodusSingleConn