function shapeFromSpecularFlow(surfChoice, pSize, KeepFlowRatio, throwBadFlow, throwParabolicCurvPoints, RegularizeSurface, verb, saveResults)

load(sprintf('../data/Surf%d_Rot10_PointCorres.mat', surfChoice));

if (throwParabolicCurvPoints)
    minLE = interp2(xGrid, yGrid, minE, X_fl, Y_fl, 'bilinear');
    maxLE = interp2(xGrid, yGrid, maxE, X_fl, Y_fl, 'bilinear');

    prPts = find(log10(maxLE./minLE) > 1.5);
    flow_x(prPts) = 0/0;
    flow_y(prPts) = 0/0;
    imagesc([flow_x flow_y]); pause(3);
end

if (throwBadFlow)
    indx = find(mScore > 1);
    flow_x(indx) = 0/0;
    flow_y(indx) = 0/0;
end

if (KeepFlowRatio ~= 1)
   len0 = length(find(isnan(flow_x)));
    
    
   len = prod(size(flow_x));
   randIndx = randperm(len - len0);
   baseIndx = setdiff(1:len, find(isnan(flow_x)));
   randIndx = baseIndx(randIndx);
   if (floor((1-KeepFlowRatio)*len) > len0)
       randIndx = randIndx(1:floor((1-KeepFlowRatio)*len-len0));
       flow_x(randIndx) = 0/0;
       flow_y(randIndx) = 0/0;
   end
end


pSizeR = pSize; pSizeC = pSizeR;
pShiftR = pSizeR; pShiftC = pSizeR;


[OvMat1, OvMat2, pStack, SMap, DMap, pMap] = generateSpecularFlowEquations_WeakNeighbors(X_fl, Y_fl, flow_x, flow_y, R, pSizeR, pSizeC, pShiftR, pShiftC, verb);
pCount = length(pStack);

OvMat3 = [];
if (RegularizeSurface)

    zIndx = find(pMap(:) <= 0    );
    OvMat3 = [];
    for ijk = 1:length(zIndx)
        tmp = zeros(2, pCount*5);
        tmp(:, (zIndx(ijk)-1)*5+[3 4]) = eye(2); %*1/(1+sqrt(pMap(ijk))*.1);
        OvMat3 = [ OvMat3; tmp ];
    end
    OvMat3 = OvMat3/10;
    OvMat3 = sparse(OvMat3);
end


OvMat = [ OvMat1; OvMat2; OvMat3 ];
OvMat_Sq = full( OvMat'*OvMat);

disp(sprintf('EVP: Size of Problem %d x %d\n', size(OvMat_Sq,1), size(OvMat_Sq, 2)));
eigSolve = 1;

if (size(OvMat_Sq, 1)) > 4000
    disp('ARE U SURE U WANT TO PROCEED !!! O_O !! Press any key to proceed. or Cntrl+c to end');
    k = input('Enter 1 if u want CVX Solver');
    if (k==1)
        eigSolve = 0;
    end
end

if (eigSolve)
    [A, B] = eig(OvMat_Sq);
    B = log10(diag(B));
    disp(B(1:10));

    [mVal, mIndx] = min(B);
    xAns = A(:, mIndx);

else
    cvx_precision best
    cvx_begin
    variables x(pCount*5)
    minimize norm( OvMat*x, 2)
    subject to
    sum(x) >= 1;
    cvx_end
    xAns = x/norm(x);
end


pPresCount = zeros(size(X_fl));
fx_recons = zeros(size(X_fl));
fy_recons = zeros(size(X_fl));
for nP = 1:pCount
    zVec = xAns((nP-1)*5+(1:5));

    tmp_fx = (zVec(1) + zVec(3)*pStack(nP).xGrid+zVec(4)*pStack(nP).yGrid);
    tmp_fy = (zVec(2) + zVec(4)*pStack(nP).xGrid+zVec(5)*pStack(nP).yGrid);

    fx_recons(pStack(nP).rIndx+(1:size(pStack(nP).xGrid,1)), pStack(nP).cIndx+((1:size(pStack(nP).xGrid,2)))) = fx_recons(pStack(nP).rIndx+(1:size(pStack(nP).xGrid,1)), pStack(nP).cIndx+((1:size(pStack(nP).xGrid,2)))) + tmp_fx;
    fy_recons(pStack(nP).rIndx+(1:size(pStack(nP).xGrid,1)), pStack(nP).cIndx+((1:size(pStack(nP).xGrid,2)))) = fy_recons(pStack(nP).rIndx+(1:size(pStack(nP).xGrid,1)), pStack(nP).cIndx+((1:size(pStack(nP).xGrid,2))))+ tmp_fy;

    pPresCount(pStack(nP).rIndx+(1:size(pStack(nP).xGrid,1)), pStack(nP).cIndx+((1:size(pStack(nP).xGrid,2)))) = pPresCount(pStack(nP).rIndx+(1:size(pStack(nP).xGrid,1)), pStack(nP).cIndx+((1:size(pStack(nP).xGrid,2)))) + 1;
end



fx_recons = fx_recons./pPresCount;
fy_recons = fy_recons./pPresCount;
scale = sum(fx_fl(:).*fx_recons(:) + fy_fl(:).*fy_recons(:))/sum(fx_recons(:).^2 + fy_recons(:).^2);
fx_recons = scale*fx_recons;
fy_recons = scale*fy_recons;


zImg1 = []; zImg2 = [];
zImg1(:, :, 1) = fx_fl; zImg1(:,:,2) = fy_fl; zImg1(:,:,3) = 1; zImg1 = zImg1./repmat( sqrt( sum(zImg1.^2, 3)), [ 1 1 3]);
zImg2(:, :, 1) = fx_recons; zImg2(:,:,2) = fy_recons; zImg2(:,:,3) = 1; zImg2 = zImg2./repmat( sqrt( sum(zImg2.^2, 3)), [ 1 1 3]);

figure
subplot 121
imshow([zImg1 zImg2]); title('Normal Fields [Original | Estimated]');
subplot 122
imagesc( real(acos(sum(zImg1.*zImg2, 3)))*180/pi); colorbar
title('Angular Error in normals in degrees');

Wflag = 0;
bimage = []; %zeros(size(Z));
TH_CURL = 1;
Z = integrability_2d(fx_recons, fy_recons, Wflag,TH_CURL,bimage);
scale = sum((f_fl(:)-mean(f_fl(:))).*(Z(:)-mean(Z(:))))/(sum((Z(:)-mean(Z(:))).^2));
Zprime = scale*(Z - mean(Z(:))) + mean(f_fl(:));

figure
subplot 121
imagesc([f_fl Zprime]); colorbar
title('Surface [Actual | Estimates ]');
subplot 122
imagesc(f_fl-Zprime); colorbar
title(sprintf('Err. Reconstruction SNR: %f', -10*log10(var(f_fl(:)-Zprime(:))/var(f_fl(:)))));
figure
mesh(f_fl); hold on; mesh(Zprime+.2);


if (saveResults)
    save Results.mat fx_recons fy_recons Zprime pSize OvMat1 OvMat2 OvMat fx_fl fy_fl f_fl xGrid yGrid f_depth f_grad_x f_grad_y zImg1 zImg2
end

keyboard
