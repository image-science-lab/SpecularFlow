function img = renderSurface_equiArea(f_depth,  xGrid, yGrid, envMap, fx, fy)

%envMap = envMap(1:size(envMap,1), 1:size(envMap,1), :);

xRes = xGrid(1,2)-xGrid(1,1);
yRes = yGrid(2,1) - yGrid(1,1);

[enX, enY] = meshgrid(1:size(envMap,2), 1:size(envMap,1));

enX  = 4*(enX/size(envMap,2))-2;
enY  = 4*(enY/size(envMap,1))-2;

if nargin < 5
    [fx, fy] = gradient(f_depth, xRes, yRes);
end

angX = -2*fx./sqrt(1+fx.^2+fy.^2);
angY = -2*fy./sqrt(1+fx.^2+fy.^2);

img = zeros([size(xGrid) size(envMap,3)]);

for ch=1:size(envMap,3)
   img(:,:, ch) = interp2(enX, enY, envMap(:,:,ch), angX, angY, 'bicubic'); 
end