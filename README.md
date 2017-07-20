# RetinalVesselSegmentationMaxJ

Authors:
Emanuele Del Sozzo (emanuele.delsozzo@polimi.it)
Marcello Pogliani (marcello.pogliani@polimi.it)

This repository contains a hardware acceleration of a retinal vessel segmentation application.
The application was accelerated on Maxeler DataFlow Engine (DFE).

This implementation won the 2nd prize at Maxeler Open Dataflow Design Competition held during the 25th IEEE International Symposium on Field-Programmable Custom Computing Machines (FCCM 2017), from April 30 to May 2, in Napa Valley, California, USA.

This repository is organized as follows:

- folder "bin" contains the executable, a README for its usage, and a test image ("original.png")
- folder "code" contains the CpuCode and EngineCode
- "report.pdf" contains the description of the work done

If you want to replicate the results of this work, please remember to add OpenMP and OpenCV flags to "Makefile.rules" file.
In particular, add: 
- "-fopenmp" to CFLAGS
- "-lcv -lhighgui" to LDFLAGS

This work has been developed and tested on a MaxWorkstation: VECTIS DFE (Max3 Card) / Intel Core i7 870 @ 2.93 GHz / 16 GB RAM
MaxCompiler 2015.2