clear all
close all


addpath('./functions/');


verb = 0; %%DEMO PARAMETER: Verbose
winSizeQuadric = 11; %%%%% DEMO PARAMETER
winSizePlanar = 5; %%%%% DEMO PARAMETER



quadricAnalysis(20, winSizeQuadric, verb);
planarAnalysis(20, winSizePlanar, verb);
disp('Press a key....');
pause; close all;
quadricAnalysis(21, winSizeQuadric, verb);
planarAnalysis(21, winSizePlanar, verb);
disp('Press a key....');
pause; close all;
quadricAnalysis(22, winSizeQuadric, verb);
planarAnalysis(22, winSizePlanar, verb);
disp('Press a key....');
pause; close all;
quadricAnalysis(23, winSizeQuadric, verb);
planarAnalysis(23, winSizePlanar, verb);

