function [ pStack, OvMat1, OvMat2, OvMat3] = generateSpecularFlowEquations_WeakNeighbors_SparseFlow(X_fl, Y_fl, X_sta, Y_sta, X_stp, Y_stp, R, pSizeR, pSizeC, pShiftR, pShiftC, verb)

RTrans = R';

r11= R(1,1); r12 = R(1,2); r21 =R(2,1); r22 = R(2,2);

nPRow = floor((size(X_fl, 1)-pShiftR+1)/pShiftR);
nPCol = floor((size(X_fl, 2)-pShiftC+1)/pShiftC);


%%ASSEMBLE PATCHES
disp('Assembling Patches');
pCount = 0;
minMax = [];
for nC = 1:nPCol
    if (verb)
        disp(nPCol - nC)
    end
    for nR=1:nPRow

        tmpMask = zeros(size(X_fl));

        XPatch_fl = X_fl( (nR-1)*pShiftR+(1:pSizeR), (nC-1)*pShiftC+(1:pSizeC) );
        YPatch_fl = Y_fl( (nR-1)*pShiftR+(1:pSizeR), (nC-1)*pShiftC+(1:pSizeC) );
        if (nR == nPRow)
            XPatch_fl = X_fl( ((nR-1)*pShiftR+1):end, (nC-1)*pShiftC+(1:pSizeC) );
            YPatch_fl = Y_fl( ((nR-1)*pShiftR+1):end, (nC-1)*pShiftC+(1:pSizeC) );
        end
        if (nC == nPCol)
            XPatch_fl = X_fl( (nR-1)*pShiftR+(1:pSizeR), ((nC-1)*pShiftC+1):end );
            YPatch_fl = Y_fl( (nR-1)*pShiftR+(1:pSizeR), ((nC-1)*pShiftC+1):end );
        end

        if (nC == nPCol) & (nR == nPRow)
            XPatch_fl = X_fl( ((nR-1)*pShiftR+1):end, ((nC-1)*pShiftC+1):end );
            YPatch_fl = Y_fl( ((nR-1)*pShiftR+1):end, ((nC-1)*pShiftC+1):end );
        end

        pCount = pCount + 1;



        zPatch.rIndx = (nR-1)*pShiftR;
        zPatch.cIndx = (nC-1)*pShiftC;
        zPatch.xGrid = XPatch_fl;
        zPatch.yGrid = YPatch_fl;

        minMax = [ minMax; min(XPatch_fl(:)) max(XPatch_fl(:)) min(YPatch_fl(:)) max(YPatch_fl(:)) ];
        pStack(pCount) = zPatch;

    end
end

disp('Neighborhood Constraints');
OvMat1 = []; %sparse(200000, pCount*5);
OvMat2 = []; %sparse(200000, pCount*5);

SMap = zeros(size(X_fl));
DMap = zeros(size(X_fl));

xFlRes = (X_fl(1,2)-X_fl(1,1));
yFlRes = (Y_fl(2,1)-Y_fl(1,1));


for row=1:size(X_fl, 1)
    if (verb)
        disp(size(X_fl,1)-row)
    end
    for col=1:size(X_fl, 2)
        %Get Source Indices
        xs = X_fl(row, col); ys = Y_fl(row, col);
        sIndx = find( ( xs >= minMax(:,1)- xFlRes ) & (xs <= minMax(:,2)+xFlRes) & (ys >= minMax(:, 3)-yFlRes) & (ys <= minMax(:,4)+yFlRes));
        SMap(row, col) = length(sIndx);
        if (isempty(sIndx))
            disp('STUPID ERROR!!!!!!');
        end

        S1Mat = [1 0 xs ys 0; 0 1 0 xs ys];

%         if (length(sIndx) > 1)
%             tmp = zeros(2*(length(sIndx)-1), 5*pCount);
%             count = 0;
%             for k1=1:length(sIndx)-1
%                 tmp(2*count+(1:2), (sIndx(k1)-1)*5+(1:5)) = S1Mat;
%                 tmp(2*count+(1:2), (sIndx(k1+1)-1)*5+(1:5)) = -S1Mat;
%                 count = count+1;
%             end
%             OvMat1 = [ OvMat1; sparse(tmp)];
% 
%         end


        xsN =  X_fl(row, col) + xFlRes/2;
        ysN = Y_fl(row, col);
        SNMat = [1 0 xsN ysN 0; 0 1 0 xsN ysN];
        sIndxN = find( ( xsN >= minMax(:,1)-xFlRes/2) & (xsN <= minMax(:,2)+xFlRes/2) & (ysN >= minMax(:, 3)-yFlRes/2) & (ysN <= minMax(:,4)+yFlRes/2));
        SMap(row, col) =  SMap(row, col) + length(sIndxN);
        if (length(sIndxN) > 1)
            tmp = zeros(2*(length(sIndxN)-1), 5*pCount);
            count = 0;
            for k1=1:length(sIndxN)-1
                tmp(2*count+(1:2), (sIndxN(k1)-1)*5+(1:5)) = SNMat;
                tmp(2*count+(1:2), (sIndxN(k1+1)-1)*5+(1:5)) = -SNMat;
                count = count+1;
            end
            OvMat1 = [ OvMat1; sparse(tmp)];

        end
        xsN =  X_fl(row, col);
        ysN = Y_fl(row, col)+ yFlRes/2;
        SNMat = [1 0 xsN ysN 0; 0 1 0 xsN ysN];
        sIndxN = find( ( xsN >= minMax(:,1)-xFlRes/2) & (xsN <= minMax(:,2)+xFlRes/2) & (ysN >= minMax(:, 3)-yFlRes/2) & (ysN <= minMax(:,4)+yFlRes/2));
        SMap(row, col) =  SMap(row, col) + length(sIndxN);
        if (length(sIndxN) > 1)
            tmp = zeros(2*(length(sIndxN)-1), 5*pCount);
            count = 0;
            for k1=1:length(sIndxN)-1
                tmp(2*count+(1:2), (sIndxN(k1)-1)*5+(1:5)) = SNMat;
                tmp(2*count+(1:2), (sIndxN(k1+1)-1)*5+(1:5)) = -SNMat;
                count = count+1;
            end
            OvMat1 = [ OvMat1; sparse(tmp)];

        end
        
        
        
    end
end

disp('Flow Constraints');
dCount = 0;
pMap = zeros(nPRow, nPCol);
for ijk=1:length(X_sta)
    if (verb)
        disp(length(X_sta)-ijk)
    end
    xs = X_sta(ijk); ys = Y_sta(ijk);
    sIndx = find( ( xs >= minMax(:,1)-xFlRes ) & (xs <= minMax(:,2)+xFlRes) & (ys >= minMax(:, 3)-yFlRes) & (ys <= minMax(:,4)+yFlRes));
    S1Mat = [1 0 xs ys 0; 0 1 0 xs ys];
    
    %Get destination Indices
    xd = X_stp(ijk); yd = Y_stp(ijk);

    xdr = RTrans(1,1)*xd + RTrans(1,2)*yd;
    ydr = RTrans(2,1)*xd + RTrans(2,2)*yd;

    dIndx = find( ( xdr >= minMax(:,1)-xFlRes) & (xdr <= minMax(:,2)+xFlRes) & (ydr >= minMax(:, 3)-yFlRes) & (ydr <= minMax(:,4)+yFlRes));

    D1Mat = [ r11 r12 r11*xdr (r11*ydr+r12*xdr) (r12*ydr);  r21 r22 r21*xdr (r21*ydr+r22*xdr) (r22*ydr)];
    
    if (~isempty(sIndx)) &  (~isempty(dIndx))
        if ~isempty(dIndx)
            %Compile matrix
            tmp = zeros(2*length(sIndx)*length(dIndx), pCount*5);
            count  = 0;
            for k1=1:length(sIndx)
                pMap(sIndx(k1)) = pMap(sIndx(k1))+1;
                for k2=1:length(dIndx)
                    tmp(2*count+(1:2), (sIndx(k1)-1)*5+(1:5)) = S1Mat;
                    tmp(2*count+(1:2), (dIndx(k2)-1)*5+(1:5)) = tmp(2*count+(1:2), (dIndx(k2)-1)*5+(1:5))-D1Mat;
                    count = count + 1;
                    
                    
                end
            end

            OvMat2 = [ OvMat2; sparse(tmp)];
            dCount = dCount + 1;
        end
    end
end

%%Prune Answer Space
%[pIR, pIC] = ind2sub(size(pMap), find(pMap(:)));
%[X, Y] = meshgrid(1:nPCol, 1:nPRow);
%pTrueMap = griddata(pIC, pIR, pMap(find(pMap(:))), X, Y, 'linear');
%pTrueMap(isnan(pTrueMap)) = 0;
%zIndx = find((pTrueMap(:) <= 2)); zIndx = 5*(zIndx(:)-1)*ones(1,5) + ones(length(zIndx(:)),1)*[1:5]; zIndx = zIndx';  zIndx = zIndx(:);
%zIndx1 = find( sum(abs(OvMat1(:,  zIndx)), 2) );
%zIndx2 = find( sum(abs(OvMat2(:,  zIndx)), 2) );

%Ov1 = OvMat1; Ov2 = OvMat2;
%Ov1(zIndx1, :) = []; Ov1(:, zIndx) = [];
%Ov2(zIndx2, :) = []; Ov2(:, zIndx) = [];
zIndx = find(pMap(:) <= 0.05*size(OvMat2,1)/prod(size(pMap)));
%zIndx= 1:length(pMap(:));
OvMat3 = [];
for ijk = 1:length(zIndx)
    tmp = zeros(3, pCount*5);
    tmp(:, (zIndx(ijk)-1)*5+(3:5)) = eye(3); %*1/(1+sqrt(pMap(ijk))*.1);
    OvMat3 = [ OvMat3; tmp ];
end
OvMat3 = OvMat3/10;
OvMat3 = sparse(OvMat3);