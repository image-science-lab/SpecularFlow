function [OvMat1, OvMat2, pStack, SMap, DMap, pMap] = generateSpecularFlowEquations_WeakNeighbors(X_fl, Y_fl, flow_x, flow_y, R, pSizeR, pSizeC, pShiftR, pShiftC, verb)

RTrans = R';

r11= R(1,1); r12 = R(1,2); r21 =R(2,1); r22 = R(2,2);

%nPRow = ceil((size(X_fl, 1)-pSizeR)/pShiftR);
%nPCol = ceil((size(X_fl, 2)-pSizeC)/pShiftC);
nPRow = floor((size(X_fl, 1))/pShiftR);
nPCol = floor((size(X_fl, 2))/pShiftC);


%%ASSEMBLE PATCHES
pCount = 0;
minMax = [];



for nC = 1:nPCol
    for nR=1:nPRow


        XPatch_fl = X_fl( (nR-1)*pShiftR + (1:pSizeR), (nC-1)*pShiftC+(1:pSizeC) );
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



disp('Pseud Overlap Mat');
OvMat1 = []; %sparse(200000, pCount*5);
OvMat2 = []; %sparse(200000, pCount*5);

SMap = zeros(size(flow_x));
DMap = zeros(size(flow_y));

xFlRes = (X_fl(1,2)-X_fl(1,1));
yFlRes = (Y_fl(2,1)-Y_fl(1,1));

pMap = zeros(nPRow, nPCol);

for row=1:size(flow_x, 1)
    if verb
        disp(size(flow_x,1)-row)
    end
    for col=1:size(flow_x, 2)
        %Get Source Indices
        xs = X_fl(row, col); ys = Y_fl(row, col);
        sIndx = find( ( xs >= minMax(:,1) ) & (xs <= minMax(:,2)) & (ys >= minMax(:, 3)) & (ys <= minMax(:,4)));


        S1Mat = [1 0 xs ys 0; 0 1 0 xs ys];

        if (length(sIndx) > 1)
            printf('ERROR');
        end
        
        
        
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
        
        if ~isnan(flow_x(row, col))
            %Get destination Indices
            xd = X_fl(row, col)+flow_x(row,col); yd = Y_fl(row, col)+flow_y(row, col);

            xdr = RTrans(1,1)*xd + RTrans(1,2)*yd;
            ydr = RTrans(2,1)*xd + RTrans(2,2)*yd;

            dIndx = find( ( xdr >= minMax(:,1)-xFlRes/2) & (xdr <= minMax(:,2)+xFlRes/2) & (ydr >= minMax(:, 3)-yFlRes/2) & (ydr <= minMax(:,4)+yFlRes/2));

            DMap(row, col) = length(dIndx);

            D1Mat = [ r11 r12 r11*xdr (r11*ydr+r12*xdr) (r12*ydr);  r21 r22 r21*xdr (r21*ydr+r22*xdr) (r22*ydr)];

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
            end
        end
    end
end