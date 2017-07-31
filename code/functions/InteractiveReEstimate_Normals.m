mFlag = 1;
xRes = xGrid(1,2)-xGrid(1,1);
yRes = yGrid(2,1)-yGrid(1,1);
[f_grad_x, f_grad_y] = gradient(f_depth, xRes, yRes);
eex = []; eey = []; ZMat = []; ZBMat = [];

minMax = [];
pCount = length(pStack);
for i = 1:length(pStack)
    minMax = [ minMax; min(pStack(i).xGrid(:)) max(pStack(i).xGrid(:)) min(pStack(i).yGrid(:)) max(pStack(i).yGrid(:)) ];
end

%OvMatN = spdiags(1./sqrt(sum(OvMat.^2, 2)), 1:size(OvMat, 1), size(OvMat,1), size(OvMat,1))*OvMat;
xFlRes = X_fl(1,2)-X_fl(1,1);
yFlRes = Y_fl(2,1)-Y_fl(1,1);


while (mFlag)
    figure(1)
    hold off
    subplot 231; imagesc(f_grad_x);
    subplot 234; imagesc(f_grad_y);
    subplot 232; imagesc(fx_recons);
    subplot 235; imagesc(fy_recons);
    subplot 233; imagesc(eex);
    subplot 236; imagesc(eey);

    subplot 231
    a = ginput;
    mFlag = 0;
    if ~isempty(a)
        mFlag = 1;
        xkno = interp2(xGrid, a(:,1), a(:,2), 'bilinear');
        ykno = interp2(yGrid, a(:,1), a(:,2), 'bilinear');
        fxno = interp2(f_grad_x, a(:,1), a(:,2), 'bilinear');
        fyno = interp2(f_grad_y, a(:,1), a(:,2), 'bilinear');
        eex = griddata(xkno, ykno, fxno, X_fl, Y_fl, 'cubic');
        eey = griddata(xkno, ykno, fyno, X_fl, Y_fl, 'cubic');
        ZMat = []; ZBMat = [];
        for ijk=1:length(xkno)
            zqIndx = find( (xkno(ijk) >= minMax(:,1)-xFlRes) & (xkno(ijk) <= minMax(:,2)+xFlRes) & (ykno(ijk) >= minMax(:,3)-yFlRes) & (ykno(ijk) <= minMax(:,4)+yFlRes) );
            if ~isempty(zqIndx)
                zqIndx = zqIndx(1);
                tmpQMat = [1 0 xkno(ijk) ykno(ijk) 0; 0 1 0 xkno(ijk) ykno(ijk) ];
                tmpBMat = [ fxno(ijk); fyno(ijk) ];
                tmpQ1Mat = zeros(2, pCount*5);
                tmpQ1Mat(:, (zqIndx-1)*5+(1:5)) = tmpQMat;
                ZMat = [ZMat; tmpQ1Mat];
                ZBMat = [ZBMat; tmpBMat];
            end
        end
        
        OvMat=[OvMat1;OvMat2];
        cvx_precision best
        cvx_begin
        variables x(pCount*5)
        minimize norm( OvMat*x, 2)
        %minimize sum(huber(OvMat*x,3))
        %minimize norm( OvMat*x, 1)
        subject to
        ZMat*x == ZBMat;
        cvx_end
        xAns = x/norm(x);

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


    end
end
