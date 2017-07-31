function planarAnalysis(surf, winSize, verb)
%%%%%%%%%%%Analytical surface descriptions
%%%SURFACE CREATION
set(0,'defaultaxesfontsize',12);

winSize = 2*floor(winSize/2)+1;
stp2 = winSize;
stp1 = floor(winSize/2);

[xyprofile, xGrid, yGrid, cutoff, analytical, f_depth ] = surfaceCreator(surf); 
if (analytical)
    x = xGrid; y = yGrid;
    f_depth = eval(xyprofile);
    if ~isempty(cutoff)
        f_depth = (f_depth >= cutoff).*(f_depth - cutoff);
    end
end
    

xRes = xGrid(1,2)-xGrid(1,1);
yRes = yGrid(2,1)-yGrid(1,1);

[f_grad_x, f_grad_y] = gradient(f_depth, xRes, yRes);


X_fl = xGrid(1+stp1:stp2:end-stp1, 1+stp1:stp2:end-stp1);
Y_fl = yGrid(1+stp1:stp2:end-stp1, 1+stp1:stp2:end-stp1);

[zX, zY] = meshgrid(-stp1:stp1, -stp1:stp1);
zX = zX*xRes;
zY = zY*yRes;

f_stack = [];
fx_stack = [];
fy_stack = [];


for row = 1:size(X_fl, 1)
    frow_stk = [];
    fxrow_stk = [];
    fyrow_stk = [];
    if (verb)
        disp(size(X_fl,1)-row)
    end
    
    for col=1:size(X_fl, 2)
        
        f_fl = f_depth( (row-1)*stp2+(1:stp2), (col-1)*stp2+(1:stp2));
        fx_fl = f_grad_x( (row-1)*stp2+(1:stp2), (col-1)*stp2+(1:stp2));
        fy_fl = f_grad_y( (row-1)*stp2+(1:stp2), (col-1)*stp2+(1:stp2));
        
        J(1) = mean(fx_fl(:));
        J(2) = mean(fy_fl(:));
        tmpfx = J(1)+0*zX;
        tmpfy = J(2)+0*zX;
        
        tmpf = J(1)*zX + J(2)*zY;
        shft = mean(f_fl(:)) - mean(tmpf(:));
        tmpf = tmpf + shft;
        
        frow_stk = [frow_stk tmpf];
        fxrow_stk = [fxrow_stk tmpfx];
        fyrow_stk = [fyrow_stk tmpfy];
        
       
    end
    f_stack = [ f_stack; frow_stk];
    fx_stack = [ fx_stack; fxrow_stk];
    fy_stack = [ fy_stack; fyrow_stk];
    
end

f_depth = f_depth(1:size(f_stack,1), 1:size(f_stack, 2));
f_grad_x = f_grad_x(1:size(f_stack,1), 1:size(f_stack, 2));
f_grad_y = f_grad_y(1:size(f_stack,1), 1:size(f_stack, 2));
xGrid = xGrid(1:size(f_stack,1), 1:size(f_stack, 2));
yGrid = yGrid(1:size(f_stack,1), 1:size(f_stack, 2));


figure
subplot 121
imagesc([f_depth f_stack]); title(['Depth Map [Ground Truth | PLANAR Approximation at block size: ' sprintf('%d', stp2) ']']);
axis equal; axis tight
subplot 122
imagesc(log10(abs((f_depth-f_stack))./(abs(f_depth)+eps))); colorbar; title('Approx. Error Ratio in log scale');
axis equal; axis tight
set(gcf,'units','normalized','outerposition',[0 0 1 1]); drawnow


zImg1 = []; zImg2 = [];
zImg1(:,:,1) = f_grad_x; zImg1(:,:,2) = f_grad_y; zImg1(:, :, 3) = 1; zImg1 = zImg1./repmat(   sqrt(sum(zImg1.^2, 3)),[1 1 3]);
zImg2(:,:,1) = fx_stack; zImg2(:,:,2) = fy_stack; zImg2(:, :, 3) = 1; zImg2 = zImg2./repmat(   sqrt(sum(zImg2.^2, 3)),[1 1 3]);
figure
subplot 121
imshow(0.5*[zImg1 zImg2]+0.5); title(['Normal Field [Ground Truth | PLANAR Approximation at block size: ' sprintf('%d', stp2) ']']);
axis equal; axis tight
subplot 122
imagesc(real(acos(sum(zImg1.*zImg2, 3)))*180/pi); colorbar
title('Angular Error in degrees');
axis equal; axis tight
set(gcf,'units','normalized','outerposition',[0 0 1 1]); drawnow

%Environmental map
load ..\data\grace-new.mat
envMap = real(hdrImg.^(1/2.2));

imgO = renderSurface_equiArea(f_depth,  xGrid, yGrid, envMap);
imgQ = renderSurface_equiArea(f_stack,  xGrid, yGrid, envMap, fx_stack, fy_stack);
figure
imshow([imgO imgQ]);
title(['Rendered Surface [Ground Truth | PLANAR Approximation at block size: ' sprintf('%d', stp2) ']']);
set(gcf,'units','normalized','outerposition',[0 0 1 1]); drawnow
