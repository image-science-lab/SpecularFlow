function img = renderSurface_euler(f_depth,  xGrid, yGrid, envMap)

xRes = xGrid(1,2)-xGrid(1,1);
yRes = yGrid(2,1) - yGrid(1,1);

envNew = zeros(size(envMap).*[1 2 1]);
envNew(:,1:size(envMap,2),:) = envMap;
envNew(:,size(envMap,2)+(1:size(envMap,2)),:) = envMap(:,end:-1:1,:);
envMap = imresize(envNew, 0.5, 'bilinear');

[enPhi, enTh] = meshgrid(1:size(envMap,2), 1:size(envMap,1));

enPhi  = 2*pi*((enPhi-1)/(size(envMap,2)-1))-pi;
enTh  = pi*((enTh-1)/size(envMap,1));


[fx, fy] = gradient(f_depth, xRes, yRes);

angTh = angle( 1 - fx.^2-fy.^2 + j*2*sqrt(fx.^2+fy.^2));
angPhi = angle(fx + j*fy);

img = zeros([size(xGrid) size(envMap,3)]);

for ch=1:size(envMap,3)
   img(:,:, ch) = interp2(enPhi, enTh, envMap(:,:,ch), angPhi, angTh, 'bicubic'); 
end