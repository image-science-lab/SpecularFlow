function [xyProfile, xGrid, yGrid, cutoff, analytical, f_depth] = surfaceCreator(choice)
cutoff = [];
analytical = 1; 
switch choice

    
    case 20, %%Simplest surface for experiments
        H = [2 -1/2; -1/2 -3/2];
        xyProfile = sprintf('(2 - (%f*(x).^2 + 2*%f.*(x).*(y)+%f*(y).^2))', H(1,1), H(1,2), H(2,2));
        [xGrid, yGrid] = meshgrid(-0.5:0.001:0.5, -0.5:0.001:0.5);
        cutoff =  [];
        
    case 21, %%Simple surface
        xyProfile = '(exp(-(1/0.52).^2*((x).^2 +2*0.67*1.414*(x).*(y)+ 2*(y).^2)))';
        [xGrid, yGrid] = meshgrid(-0.5:0.001:0.5, -0.5:0.001:0.5);
        cutoff = [];
    case 22, %%MORE COMPLEX LOOKING STUFF
        xyProfile = '0.85*(exp(- (1/0.35)^2*((x+0.3).^2 +2*0.67*1.414*(x+0.3).*(y+0.1)+ 2*(y+0.1).^2)) + 0.9*exp(- (1/0.65)^2*(4*(x-0.3).^2 + (y+0.1).^2)))';
       [xGrid, yGrid] = meshgrid(-0.5:0.001:.5, -0.5:0.001:0.5);
       cutoff = [];
    case 23,
        %%%THIS IS A REAL SURFACE - BUDDHA
        xyProfile=[];
        load ../data/Surf23_Rot10_PointCorres f_depth xGrid yGrid
        cutoff = [];
        analytical = 0;
    case 41, %%MORE COMPLEX LOOKING STUFF
        xyProfile = '(1/3)*sqrt(max(0, 4-(x*2).^2-(y*2).^2-cos(2*x*2-2)-sin(2*y*2)))';
       [xGrid, yGrid] = meshgrid(-0.5:0.001:.5, -0.5:0.001:0.5);
       cutoff = [];
    case 42, %%MORE COMPLEX LOOKING STUFF
        xyProfile = '(1/3)*sqrt(max(0, 4-(x*2).^2-(y*2).^2-cos(3*x*2-6)-2*sin(2*y*2)))';
       [xGrid, yGrid] = meshgrid(-0.5:0.001:.5, -0.5:0.001:0.5);
       cutoff = [];
    
    case 3,
        xProfile = '((abs(x) >= 0.8).*0 + (abs(x) < 0.8).*( (x+0.8).*sqrt(0.1*abs(-x+0.8))))';
        yProfile = '((abs(y) >= 0.8).*0 + (abs(y) < 0.8).*sqrt(- y.^2 + (0.8.^2) ))';   
        xyProfile = [ xProfile '.*' yProfile ];
        [xGrid, yGrid] = meshgrid(-1:0.0025:1, -1:0.0025:1);
        
    case 4,
        %%REASONABLY SINPLE ... HAS PARABOLIC POINT ISSUES
        xProfile = '( 0.75*exp(-((x+0.3)/0.3).^2)+ exp(-((x-0.2)/0.35).^2))';
        yProfile = '((abs(y) >= 0.8).*0 + (abs(y) < 0.8).*(- y.^2 + (0.8.^2) ))';
        xyProfile = [ xProfile '.*' yProfile ];
        [xGrid, yGrid] = meshgrid(-1:0.0025:1, -1:0.0025:1);
        
    case 1,
        %% NO parabolic --- Quadric
        H = [2 -1/2; -1/2 3/2];
        xyProfile = sprintf('(2 - (%f*(x).^2 + 2*%f.*(x).*(y)+%f*(y).^2))', H(1,1), H(1,2), H(2,2));
        [xGrid, yGrid] = meshgrid(-1.2:0.002:1.2, -1.2:0.002:1.2);
        cutoff =  0.5;
    
    case 2, %%SIMPLESE PARABOLIC --- GAUSSIAN
        xyProfile = '(exp(-(1/0.75).^2*((x).^2 +2*0.67*1.414*(x).*((y)+0.1)+ 2*((y)+0.1).^2)))';
        [xGrid, yGrid] = meshgrid(-1.2:0.002:1.2, -1.2:0.002:1.2);
        cutoff = 0.33;

    case 5, %%MORE COMPLEX LOOKING STUFF
        xyProfile = '0.85*(exp(- (1/0.35)^2*((x+0.3).^2 +2*0.67*1.414*(x+0.3).*(y+0.1)+ 2*(y+0.1).^2)) + 0.9*exp(- (1/0.65)^2*(4*(x-0.3).^2 + (y+0.1).^2)))';
       [xGrid, yGrid] = meshgrid(-1.2:0.002:1.2, -1.2:0.002:1.2);
       cutoff = 0.13;
       
    case 6, %%MORE COMPLEX LOOKING STUFF
        xyProfile = '0.85*(exp(- (1/0.35)^2*((x+0.3).^2 +2*0.67*1.414*(x+0.3).*(y+0.1)+ 2*(y+0.1).^2)) + 0.9*exp(- (1/0.65)^2*(4*(x-0.3).^2 + (y+0.1).^2)))';
       [xGrid, yGrid] = meshgrid(0:0.002:0.8, -0.6:0.002:0.6);
       cutoff = [];
    case 7, %%MORE COMPLEX LOOKING STUFF
        xyProfile = '0.85*(exp(- (1/0.35)^2*((x+0.3).^2 +2*0.67*1.414*(x+0.3).*(y+0.1)+ 2*(y+0.1).^2)) + 0.9*exp(- (1/0.65)^2*(4*(x-0.3).^2 + (y+0.1).^2)))';
       [xGrid, yGrid] = meshgrid(-0.5:0.002:.4, -0.5:0.002:0.2);
       cutoff = [];
       
    case 8, %%SIMPLESE PARABOLIC --- GAUSSIAN
        xyProfile = '(exp(-(1/0.75).^2*((x).^2 +2*0.67*1.414*(x).*(y+0.1)+ 2*(y+0.1).^2)))';
        [xGrid, yGrid] = meshgrid(-1.2:0.004:1.2, -1.2:0.004:1.2);
%         [xGrid, yGrid] = meshgrid(-0.6:0.004:0.6, -0.6:0.004:0.6);
        cutoff = [];
    case 9, %%SIMPLESE PARABOLIC --- GAUSSIAN
        xyProfile = '(exp(-(1/0.75).^2*((x).^2 +2*0.67*1.414*(x).*(y+0.1)+ 2*(y+0.1).^2)))';
%        [xGrid, yGrid] = meshgrid(-1.2:0.004:1.2, -1.2:0.004:1.2);
%         [xGrid, yGrid] = meshgrid(-0.6:0.004:0.6, -0.6:0.004:0.6);
        [xGrid, yGrid] = meshgrid(-1.2:0.001:1.2, -1.2:0.001:1.2);
%        cutoff = [];

    case 10,
%        xyProfile = 'exp(-(sin(5*x)).^2-(cos(5*y)).^2)';
         xyProfile = ' sin((x.^2+y.^2)*20).*exp(x)';
%         xyProfile = '(exp(-25*(x+0.2).^2)+exp(-50*(x-.25).^2)).*(exp(-10*(y+0.3).^2) - exp(-40*(y-0.25).^2))';
         [xGrid, yGrid] = meshgrid(-0.5:0.001:0.5, -0.5:0.001:0.5);
         cutoff = [];
    case 11,
         xyProfile = 'x.^3 + y.^2';
        [xGrid, yGrid] = meshgrid(-0.5:0.001:0.5, -0.5:0.001:0.5);
        cutoff =  [];
        
    case 12, %%Simplest surface for experiments
        H = [2 -1/2; -1/2 1/2];
        xyProfile = sprintf('(2 - (%f*(x).^2 + 2*%f.*(x).*(y)+%f*(y).^2))', H(1,1), H(1,2), H(2,2));
        [xGrid, yGrid] = meshgrid(-0.5:0.001:0.5, -0.5:0.001:0.5);
        cutoff =  [];
end

if (analytical)
    x= xGrid; y = yGrid;
    f_depth = eval(xyProfile);
end