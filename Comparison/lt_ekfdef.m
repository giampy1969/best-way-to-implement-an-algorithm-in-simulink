% This MATLAB script initializes a legacy code tool object in the workspace
% and then creates an S-function, MEX and TLC files, and optionally a Simulink block,
% which implements a 2-state EKF for attitude estimation
% For more info on the algorithm, its implementation and usage, go to 
% the root folder and see the files README.txt, Presentation.pdf, and 
% ekf2sim_example.slx (read the help of the ekf2sim block therein). 
% Copyright 2016 The MathWorks, Inc.

ekfdef=legacy_code('initialize');

ekfdef.SFunctionName='ekf2lt';

ekfdef.SourceFiles={fullfile(pwd,'ekf2lt_src.c')}; % note: use fullfile !!
ekfdef.SampleTime=T;

ekfdef.InitializeConditionsFcnSpec='void ekf_init(double work1[9], double work2[33], double p5[2][2], double p6[2], double p7[3])';
ekfdef.OutputFcnSpec='void ekf_out(double y1[3], double u1[3], double u2[3], double u3[3], double work1[9], double work2[33], double p1[2][2], double p2[2][2], double p3, double p4)';

legacy_code('sfcn_cmex_generate',ekfdef)
legacy_code('sfcn_tlc_generate',ekfdef)
legacy_code('compile',ekfdef)
%legacy_code('slblock_generate',ekfdef)
