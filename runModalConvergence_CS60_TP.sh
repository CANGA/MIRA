# Run self convergence test on Topography, constant mesh spacing, variable number of modes
python CANGAFieldGenDriver.py --pm meshes/outCSne60.g --so 2 --nm 32 --EvaluateTPO --ExodusSingleConn
python CANGAFieldGenDriver.py --pm meshes/outCSne60.g --so 2 --nm 64 --EvaluateTPO --ExodusSingleConn
python CANGAFieldGenDriver.py --pm meshes/outCSne60.g --so 2 --nm 128 --EvaluateTPO --ExodusSingleConn
python CANGAFieldGenDriver.py --pm meshes/outCSne60.g --so 2 --nm 256 --EvaluateTPO --ExodusSingleConn
python CANGAFieldGenDriver.py --pm meshes/outCSne60.g --so 2 --nm 512 --EvaluateTPO --ExodusSingleConn