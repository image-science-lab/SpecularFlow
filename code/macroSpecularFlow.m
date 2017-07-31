clear all
close all

addpath('./functions');
addpath('./AgrawalICCV05MatlabCode');

surfChoice = 23; %%%DEMO PARAMETER: Surface choice (20 = Quadric, 21 = Gaussian, 22 = Mixture of Gaussian, 23 = Buddha head)
KeepFlowRatio = .5; %%DEMO PARAMETER: Ratio of flow field to keep

throwParabolicCurvPoints = 0; %%DEMO PARAMETER
RegularizeSurface = 0; %%Demo Parameter
pSize = 10; %%% DEMO PARAMETER: Patch size
throwBadFlow = 1; %%DEMO PARAMETER:

verb = 0; %%% DEMO PARAMETER: Verbose
saveResults = 0; %%DEMO PARAMETER

shapeFromSpecularFlow(surfChoice, pSize, KeepFlowRatio, throwBadFlow, throwParabolicCurvPoints, RegularizeSurface, verb, saveResults);
