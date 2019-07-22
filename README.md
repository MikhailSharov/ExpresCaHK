# ExpresCaHK
Expres Ca II H and K Calibration
This project is used in order to calibrate the Calcium II H and K S values of the EXPRES instrument to those from the Mount Wilson Observatory.

The main pipeline is master_pipeline.py and should be the only file that is primarily used for calculations. Other files are supporting files, which for the most part, with a few excpetions, have been integrated into the final pipeline.
Notable exceptions are finding the RMS of fit, and any error analysis functions have not been included in the pipeline.
These must be done outside of the pipeline for now, however are not essential in understanding the meaning of the S_EXPRES.
