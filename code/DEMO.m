clear all
close all

warning off
addpath('./functions');
addpath('./AgrawalICCV05MatlabCode');

EnvMap = 5; %%%Demo PARAMETER: Environment Value %%CHOOSE DIFFERNECE VALUES GHERE FOR DIFFERENT IMAGES FROM THE "../images" folder
interactive = 0; %DEMO PARAMETER Interactive Normals. Set this to "1" 
NormalsUnif = 0; %DEMO PARAMETER Uniform no of Normals. Set to 16 
pSize = 5; %%% DEMO PARAMETER: Patch size
throwBadFlow = 1; %%DEMO PARAMETER:

verb = 0; %%% DEMO PARAMETER: Verbose
saveResults = 0; %%DEMO PARAMETER


%%%%%QUADRIC
disp('Working on QUADRIC');
surfChoice = 20; %%%DEMO PARAMETER: Surface choice (20 = Quadric)
RotationA = [ 0 ]; % 25 30 35 40 45 50]; %%Demo Parameter Rotataions
RotationShift = [ 5 ]; %%Demo Parameter Rotataions
Rotation = [RotationA(:) RotationA(:)+RotationShift]';

if (1)
    createSyntheticImages(surfChoice, EnvMap, unique(Rotation(:)));
end

close all
[X_fl, Y_fl, f_fl, fx_fl, fy_fl, Zprime, fx_recons, fy_recons]= shapeFromSpecularFlow_SIFT(surfChoice, EnvMap, Rotation, pSize, throwBadFlow, NormalsUnif, interactive, verb, saveResults);
disp('press any key to continue....');
pause
close all
drawnow

%%%%%GAUSSIAN
disp('Working on GAUSSIAN');
surfChoice = 21; %%%DEMO PARAMETER: Surface choice (21 = Gaussian)
RotationA = [ 0 5 10 15 20]; % 25 30 35 40 45 50]; %%Demo Parameter Rotataions
RotationShift = [ 5 ]; %%Demo Parameter Rotataions
Rotation = [RotationA(:) RotationA(:)+RotationShift]';

if (1)
    createSyntheticImages(surfChoice, EnvMap, unique(Rotation(:)));
end

close all
[X_fl, Y_fl, f_fl, fx_fl, fy_fl, Zprime, fx_recons, fy_recons]= shapeFromSpecularFlow_SIFT(surfChoice, EnvMap, Rotation, pSize, throwBadFlow, NormalsUnif, interactive, verb, saveResults);
disp('press any key to continue....');
pause
close all
drawnow

%%%%%complex surface 1
disp('Working on COMPLEX SURFACE 1');
surfChoice = 41; %%%DEMO PARAMETER: Surface choice (41 = Complex Surface 1)
RotationA = [ 0 5 10 ]; % 25 30 35 40 45 50]; %%Demo Parameter Rotataions
RotationShift = [ 5 ]; %%Demo Parameter Rotataions
Rotation = [RotationA(:) RotationA(:)+RotationShift]';

if (1)
    createSyntheticImages(surfChoice, EnvMap, unique(Rotation(:)));
end

close all
[X_fl, Y_fl, f_fl, fx_fl, fy_fl, Zprime, fx_recons, fy_recons]= shapeFromSpecularFlow_SIFT(surfChoice, EnvMap, Rotation, pSize, throwBadFlow, NormalsUnif, interactive, verb, saveResults);
disp('press any key to continue....');
pause
close all
drawnow



%%%%%Complex surface 2
disp('Working on COMPLEX SURFACE 2');
surfChoice = 42; %%%DEMO PARAMETER: Surface choice (41 = Complex Surface 2)
RotationA = [ 0 5 10 ]; % 25 30 35 40 45 50]; %%Demo Parameter Rotataions
RotationShift = [ 5 ]; %%Demo Parameter Rotataions
Rotation = [RotationA(:) RotationA(:)+RotationShift]';

if (1)
    createSyntheticImages(surfChoice, EnvMap, unique(Rotation(:)));
end

close all
[X_fl, Y_fl, f_fl, fx_fl, fy_fl, Zprime, fx_recons, fy_recons]= shapeFromSpecularFlow_SIFT(surfChoice, EnvMap, Rotation, pSize, throwBadFlow, NormalsUnif, interactive, verb, saveResults);
disp('press any key to continue....');
pause
close all
drawnow



%%%%%Mixture of Gaussian
disp('Working on MIXTURE OF GAUSSIAN');
disp('THIS USES ADDITIONAL NORMALS AT 16 LOCATIONS FOR BETTER ESTIMATION!!!!');
NormalsUnif = 16; %DEMO PARAMETER Uniform no of Normals. Set to 16 

surfChoice = 22; %%%DEMO PARAMETER: Surface choice (22 = MoG)
RotationA = [ 0 5 10 15 20]; % 25 30 35 40 45 50]; %%Demo Parameter Rotataions
RotationShift = [ 5 ]; %%Demo Parameter Rotataions
Rotation = [RotationA(:) RotationA(:)+RotationShift]';


if (1)
    createSyntheticImages(surfChoice, EnvMap, unique(Rotation(:)));
end

close all
[X_fl, Y_fl, f_fl, fx_fl, fy_fl, Zprime, fx_recons, fy_recons]= shapeFromSpecularFlow_SIFT(surfChoice, EnvMap, Rotation, pSize, throwBadFlow, NormalsUnif, interactive, verb, saveResults);
disp('press any key to continue....');
pause
close all
drawnow
