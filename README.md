# best-way-to-implement-an-algorithm-in-simulink
Eight ways to implement an Extended Kalman Filter as a Simulink&reg; block

[![View best-way-to-implement-an-algorithm-in-simulink on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/46786-best-way-to-implement-an-algorithm-in-simulink)

[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=giampy1969/best-way-to-implement-an-algorithm-in-simulink)

This package contains some examples and a presentation (given at the International Conference on Robotics and Automation, Hong Kong, June 2014) discussing several possible ways of implementing an algorithm in Simulink.

Specifically, a simple Extended Kalman Filter based algorithm for attitude estimation is implemented in Simulink using S-functions (in C and MATLAB), System objects:tm:, S-Function Builder, Legacy Code Tool, and the MATLAB&reg; function block (using both internal and external states).

Advantages and drawbacks of the different methods are discussed, and performance is then compared in several ways. First, the different models are simulated in Simulink, then, executable files generated from each models are executed both on an Intel laptop and on an Arduino Uno, with interesting results.
