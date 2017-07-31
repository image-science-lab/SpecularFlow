function [X_fl, Y_fl, f_fl, fx_fl, fy_fl, Zprime, fx_recons, fy_recons] =  shapeFromSpecularFlow_SIFT(surfChoice, imChoice, rotAngle, pSize, throwBadFlow, NormalsUnif, interactive, verb, saveResults)

[xyprofile, xGrid, yGrid, cutoff, analytical, f_depth] = surfaceCreator(surfChoice);
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


omega =  [0 0.01; -0.01 0]*(rotAngle(2,1) - rotAngle(1, 1));
R = expm(omega); Rt = R';

X_sta = []; Y_sta = [];
X_stp = []; Y_stp = [];
disE = [];

disp('Getting SIFT Matches'); 

for imC = imChoice
    for rotC = 1:size(rotAngle, 2)
        
        img0 = imread( sprintf('../renderedImages/Surf%02d_Rot%02d_Env%02d.png', surfChoice, rotAngle(1, rotC), imC));
        img1 = imread( sprintf('../renderedImages/Surf%02d_Rot%02d_Env%02d.png', surfChoice, rotAngle(2, rotC), imC));

        
        omega1 =  [0 0.01; -0.01 0]*rotAngle(1, rotC);
        R1 = expm(omega1); Rt1 = R1';

        if (analytical)
            x  = Rt1(1,1)*xGrid + Rt1(1,2)*yGrid;
            y  = Rt1(2,1)*xGrid + Rt1(2,2)*yGrid;
            zf_depth = eval(xyprofile);
            if ~isempty(cutoff)
                zf_depth = (zf_depth >= cutoff).*(zf_depth - cutoff);
            end

            omega1 =  [0 0.01; -0.01 0]*rotAngle(2, rotC);
            R1 = expm(omega1); Rt1 = R1';

            x  = Rt1(1,1)*xGrid + Rt1(1,2)*yGrid;
            y  = Rt1(2,1)*xGrid + Rt1(2,2)*yGrid;
            zfrot_depth = eval(xyprofile);
            if ~isempty(cutoff)
                zfrot_depth = (zfrot_depth >= cutoff).*(zfrot_depth - cutoff);
            end


            [zf_grad_x, zf_grad_y] = gradient(zf_depth, xRes, yRes);
            [zfnew_grad_x, zfnew_grad_y] = gradient(zfrot_depth, xRes, yRes);
        end
        
        

        for ch=1:3
            %Compute SIFT Features
            [F1, D1] = vl_sift(single(double(img0(:,:,ch))));
            [F2, D2] = vl_sift(single(double(img1(:,:,ch))));

            fIndx = (F1(3, :) > 10);
            F1(:, fIndx) = []; D1(:, fIndx) = [];
            fIndx = (F2(3, :) > 10);
            F2(:, fIndx) = []; D2(:, fIndx) = [];


            %Matching
            [matc, sco] = vl_ubcmatch(D1, D2, 5);

            %Displau
            figure(1)
            hold off
            imshow([img0 img1])
            hold on
            indx = random('unid', size(matc,2), 100,1); indx = unique(indx);
            plot([F1(1, matc(1,indx)); F2(1, matc(2,indx))+size(img0,2)], [F1(2, matc(1,indx)); F2(2, matc(2,indx))], 'ro-')

            drawnow

            % %Check Normals
            %
            
            n1x = -zf_grad_x./sqrt(1+zf_grad_x.^2 + zf_grad_y.^2);
            n1y = -zf_grad_y./sqrt(1+zf_grad_x.^2 + zf_grad_y.^2);
            n1z = 1./sqrt(1+zf_grad_x.^2 + zf_grad_y.^2);

            n2x = -zfnew_grad_x./sqrt(1+zfnew_grad_x.^2 + zfnew_grad_y.^2);
            n2y = -zfnew_grad_y./sqrt(1+zfnew_grad_x.^2 + zfnew_grad_y.^2);
            n2z = 1./sqrt(1+zfnew_grad_x.^2 + zfnew_grad_y.^2);

            tmp1x = interp2(n1x, F1(1, matc(1,:)), F1(2,matc(1,:)), 'bilinear');
            tmp1y = interp2(n1y, F1(1, matc(1,:)), F1(2,matc(1,:)), 'bilinear');
            tmp1z = interp2(n1z, F1(1, matc(1,:)), F1(2,matc(1,:)), 'bilinear');

            tmp2x = interp2(n2x, F2(1, matc(2,:)), F2(2,matc(2,:)), 'bilinear');
            tmp2y = interp2(n2y, F2(1, matc(2,:)), F2(2,matc(2,:)), 'bilinear');
            tmp2z = interp2(n2z, F2(1, matc(2,:)), F2(2,matc(2,:)), 'bilinear');

            dis = tmp1x.*tmp2x + tmp1y.*tmp2y +  tmp1z.*tmp2z;
            dis = acos(dis)*180/pi;
            figure(2)
            hist(dis);
            drawnow



            %Compute Flow
            X_sta1 = interp2(xGrid, F1(1, matc(1,:)), F1(2, matc(1,:)), 'bicubic'); 
            Y_sta1 = interp2(yGrid, F1(1, matc(1,:)), F1(2, matc(1,:)), 'bicubic'); 
            X_stp1 = interp2(xGrid, F2(1, matc(2,:)), F2(2, matc(2,:)), 'bicubic');
            Y_stp1 = interp2(yGrid, F2(1, matc(2,:)), F2(2, matc(2,:)), 'bicubic');
            
            %Unwarp
            omega1 =  [0 0.01; -0.01 0]*rotAngle(rotC);
            R1 = expm(omega1); Rt1 = R1';
            tmpX = Rt1*[X_sta1(:)'; Y_sta1(:)'];
            X_sta1 = tmpX(1,:); Y_sta1 = tmpX(2,:);
            tmpX = Rt1*[X_stp1(:)'; Y_stp1(:)'];
            X_stp1 = tmpX(1,:); Y_stp1 = tmpX(2,:);
            

            X_sta = [X_sta; X_sta1(:)]; Y_sta = [Y_sta; Y_sta1(:)];
            X_stp = [X_stp; X_stp1(:)];  Y_stp = [Y_stp; Y_stp1(:)];
            
            disE = [disE; dis(:)];
        end
    end
end
dis  = disE;

figure
plot(X_sta(:), Y_sta(:), '.'); axis([-0.5 0.5 -0.5 0.5]); axis ij

if (throwBadFlow)
    indx = find(dis > 1);
    X_sta(indx) = [];
    Y_sta(indx) = [];
    X_stp(indx) = [];
    Y_stp(indx) = [];
    dis(indx) = [];
end

[X_fl, Y_fl] = meshgrid(-0.5:0.02:0.5, -0.5:0.02:0.5);

pp = pSize;
pSizeR = pp; pSizeC = pp;
pShiftR = pp; pShiftC = pp;

[ pStack, OvMat1, OvMat2, OvMat3] = generateSpecularFlowEquations_WeakNeighbors_SparseFlow(X_fl, Y_fl, X_sta, Y_sta, X_stp, Y_stp, R, pSizeR, pSizeC, pShiftR, pShiftC, verb);
pCount = length(pStack);


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

f_fl = interp2(xGrid, yGrid, f_depth, X_fl, Y_fl, 'bilinear');
fx_fl = interp2(xGrid, yGrid, f_grad_x, X_fl, Y_fl, 'bilinear');
fy_fl = interp2(xGrid, yGrid, f_grad_y, X_fl, Y_fl, 'bilinear');


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



if (NormalsUnif)
    UniformlyReEstimate_Normals
end
if (interactive)
    InteractiveReEstimate_Normals
end

if (NormalsUnif) | (interactive)
    disp('USING Additional Known Normals');

    scale = sum(fx_fl(:).*fx_recons(:) + fy_fl(:).*fy_recons(:))/sum(fx_recons(:).^2 + fy_recons(:).^2);
    fx_recons = scale*fx_recons;
    fy_recons = scale*fy_recons;


    zImg1 = []; zImg2 = [];
    zImg1(:, :, 1) = fx_fl; zImg1(:,:,2) = fy_fl; zImg1(:,:,3) = 1; zImg1 = zImg1./repmat( sqrt( sum(zImg1.^2, 3)), [ 1 1 3]);
    zImg2(:, :, 1) = fx_recons; zImg2(:,:,2) = fy_recons; zImg2(:,:,3) = 1; zImg2 = zImg2./repmat( sqrt( sum(zImg2.^2, 3)), [ 1 1 3]);

    figure
    subplot 121
    imshow([zImg1 zImg2]); title('Normal Fields [Original | Estimated] (with some known normals)');
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
    title('Surface [Actual | Estimates ] (with some known normals)');
    subplot 122
    imagesc(f_fl-Zprime); colorbar
    title(sprintf('Err. Reconstruction SNR: %f', -10*log10(var(f_fl(:)-Zprime(:))/var(f_fl(:)))));
    figure
    mesh(f_fl); hold on; mesh(Zprime+.2); title('with a few known normals')
end


if (saveResults)
    save Results.mat fx_recons fy_recons Zprime pSize OvMat1 OvMat2 OvMat3 OvMat fx_fl fy_fl f_fl xGrid yGrid f_depth f_grad_x f_grad_y zImg1 zImg2 X_sta Y_sta X_stp Y_stp dis
end

