EKF: Eight ways to implement an Extended Kalman Filter as a Simulink(R) block
-----------------------------------------------------------------------------

This package contains some examples and a presentation (given at the 
International Conference on Robotics and Automation, Hong Kong, June 2014)
discussing several possible ways of implementing an algorithm in Simulink.

Specifically, a simple Extended Kalman Filter based algorithm for attitude 
estimation is implemented in Simulink using S-functions (in C and MATLAB), 
System objects(TM), S-Function Builder, Legacy Code Tool, and the 
MATLAB(R) function block (using both internal and external states).

Advantages and drawbacks of the different methods are discussed, and performance
is then compared in several ways. First, the different models are simulated
in Simulink, then, executable files generated from each models are executed
both on an Intel laptop and on an Arduino Uno, with interesting results.


Contents:
---------

README.txt              This file
license.txt             License file

Presentation.pdf	Presentation (for ICRA 2014) slides

ekf2sim_example.slx     2-states EKF example and explanation (see block help)
ekf9mlf_example.slx     9-states EKF example and explanation (see block help)
ekf9mls_example.slx     9-states EKF example and explanation (see block help)
signals.mat             File containing signals for the 2 examples above
a9_sym.m                Symbolic derivation of the A matrix for the 9-states EKF

Comparison              Folder containing the different EKF implementations
                        as discussed in the presentation.


Contents of the Comparison folder:
----------------------------------

ekf2sfm.slx             EKF as a MATLAB S-function block
sfall2.m                MATLAB S-function code for ekf2sfm

ekf2sys.slx             EKF as a MATLAB System object block
sysall.m                MATLAB code (class definition) for ekf2sys

ekf2mlf.slx             EKF as a MATLAB function block

ekf2mls.slx             EKF as a MATLAB function block with external states

ekf2sfc.slx             EKF as a C S-function block
ekf2t.c                 C S-function code for ekf2sfm

ekf2sfb.slx             EKF as an S-Function Builder block

ekf2lct.slx             EKF as a Legacy Code Tool (LCT) block
lt_ekfdef.m             MATLAB code to create the LCT block for ekf2lct
ekf2lt_src.c            C source code for ekf2lct

ekf2sim.slx             EKF assembled from basic Simulink blocks


Running the comparison models (in the Comparison folder):
---------------------------------------------------------

You must double click on the "Load Data" pink button before executing 
models in the Comparison folder. 

Also, you need to generate MEX files (only once) before running simulations in which 
the EKF is implemented with a C-based method (C S-function, S-Function Builder, LCT)

Specifically, you must:

    -) Compile the C S-function using the command "mex ekf2t.c" before 
       running ekf2sfc.slx the first time

    -) Run lt_ekfdef.m from the MATLAB command line before running
       ekf2lct.slx the first time

    -) Open the S-Function builder block and click on the "Build" button
       before running ekf2sfb.slx the first time


Versions:
---------
Version 1.0, 22-May-2014, Initial version.
Version 1.1, 28-May-2014, Added info to Simulink files and other minor refinements.
Version 1.2, 28-May-2014, Streamlined signal generator.


Trademark Notice:
-----------------
MATLAB and Simulink are registered trademarks of The MathWorks, Inc.  
See www.mathworks.com/trademarks for a list of additional trademarks.  

Copyright 2016, The MathWorks, Inc.