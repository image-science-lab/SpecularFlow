clear all
close all

addpath('./functions');
%%%%%%%%%%%Analytical surface descriptions
%%%SURFACE CREATION

[xyprofile, xGrid, yGrid, cutoff, analytical, f_depth] = surfaceCreator(23);

if (analytical)
    x = xGrid; y = yGrid;
    f_depth = eval(xyprofile);
    if ~isempty(cutoff)
        f_depth = (f_depth >= cutoff).*(f_depth - cutoff);
    end
end
imagesc(f_depth); colorbar

xRes = xGrid(1,2)-xGrid(1,1);
yRes = yGrid(2,1)-yGrid(1,1);

[f_grad_x, f_grad_y] = gradient(f_depth, xRes, yRes);


return

%%%%%%%SURFACE ANALYSIS

xFl  =xGrid; yFl = yGrid; %%Highres
%xFl  =xGrid(1:10:end, 1:10:end); yFl = yGrid(1:10:end, 1:10:end);

[ eigL, eigH] = parabolicPointAnalysis(xGrid, yGrid, f_depth, xFl, yFl, 10);
minE = min(abs(eigL), abs(eigH));
maxE = max(abs(eigL), abs(eigH));
imagesc(log10([minE maxE]),[-2 3]); colorbar
figure
imagesc(log10([maxE./minE]),[-2 4]); colorbar

return

%%%ROtate and Estimate
omega = [0 0.01; -0.01 0]*10;
R = expm(omega); Rt = R';

if (analytical)
    x = Rt(1,1)*xGrid + Rt(1,2)*yGrid;
    y = Rt(2,1)*xGrid + Rt(2,2)*yGrid;
    frot_depth = eval(xyprofile);
    if ~isempty(cutoff)
        frot_depth = (frot_depth >= cutoff).*(frot_depth - cutoff);
    end
else
    x = Rt(1,1)*xGrid + Rt(1,2)*yGrid;
    y = Rt(2,1)*xGrid + Rt(2,2)*yGrid;
    frot_depth = interp2(xGrid, yGrid, f_depth, x, y, 'bicubic', -100);
end

imagesc(frot_depth); colorbar
return


%%Estimate Specular Flow
X_fl = xGrid(1:10:end, 1:10:end);
Y_fl = yGrid(1:10:end, 1:10:end);


fx_fl = interp2(xGrid, yGrid, f_grad_x, X_fl, Y_fl, 'bilinear');
fy_fl = interp2(xGrid, yGrid, f_grad_y, X_fl, Y_fl, 'bilinear');
f_fl = interp2(xGrid, yGrid, f_depth, X_fl, Y_fl, 'bilinear');


%[flow_x, flow_y, mScore] = computeSpecularFlow(xGrid, yGrid, f_depth, frot_depth, R, X_fl, Y_fl);
[flow_x, flow_y, mScore] = computeSpecularPointMatches(xGrid, yGrid, f_depth, frot_depth, R, X_fl, Y_fl);
%[flow_x, flow_y, mScore, flow_xi, flow_yi, matchScore] = computeSpecularPointMatches_InfiniteResolution(xGrid, yGrid, f_depth, frot_depth, xyprofile, R,  X_fl, Y_fl);
%flow_x = flow_xi; flow_y = flow_yi; mScore = matchScore;
%
%

