function [eigL, eigH] = parabolicPointAnalysis(xGrid, yGrid, f_depth, xFl, yFl, nSize)

xRes = xGrid(1,2)-xGrid(1,1);
yRes = yGrid(2,1)-yGrid(1,1);

[zX, zY] = meshgrid(-nSize:nSize, -nSize:nSize);
zX = zX*xRes;
zY = zY*yRes;

[fx, fy] = gradient(f_depth, xRes, yRes);

for row=1:size(xFl, 1)
    disp(size(xFl,1)-row)
    for col=1:size(yFl, 2)


        
        idx = find(abs(xGrid(1,:) - xFl(row, col)) < 1e-5);
        idy = find(abs(yGrid(:,1) - yFl(row, col)) < 1e-5);
        
        if isempty(idx) || isempty (idy)
            xVal = xFl(row, col)+zX;
            yVal = yFl(row, col)+zY;
            fVal_x = interp2(xGrid, yGrid, fx, xVal, yVal, 'nearest');
            fVal_y = interp2(xGrid, yGrid, fy, xVal, yVal, 'nearest');
            fVal_x = fVal_x(:); xVal = xVal(:); yVal = yVal(:);        fVal_y = fVal_y(:);
            xVal(isnan(fVal_x)) = [];
            yVal(isnan(fVal_x)) = [];
            fVal_x(isnan(fVal_x)) = [];
            fVal_y(isnan(fVal_y)) = [];
        else
            indexX = max(1, idx-nSize):min(size(xGrid, 2), idx+nSize);
            indexY = max(1, idy-nSize):min(size(yGrid, 1), idy+nSize);

            xVal = xGrid(indexY, indexX);
            yVal = yGrid(indexY, indexX);
            fVal_x = fx(indexY, indexX);
            fVal_y = fy(indexY, indexX);
            fVal_x = fVal_x(:); xVal = xVal(:); yVal = yVal(:);        fVal_y = fVal_y(:);
            xVal(isnan(fVal_x)) = [];
            yVal(isnan(fVal_x)) = [];
            fVal_y(isnan(fVal_x)) = [];
            fVal_x(isnan(fVal_x)) = [];
            
        end
        
        if isnan(f_depth(idy, idx))
            eigL(row, col) = -Inf;
            eigH(row, col) = Inf;
        else

            Amat = [ xVal yVal 0*xVal 1+0*xVal 0*xVal;
                0*xVal xVal yVal 0*xVal 1+0*xVal];
            Bmat = [ fVal_x; fVal_y];

            param = Amat\Bmat;

            qParams(row, col, :) = param;

            Hmat = [ param(1) param(2); param(2) param(3)];
            [A, B] = eig(Hmat); B = diag(B);
            eigL(row, col) = min(B);
            eigH(row, col) = max(B);
        end
    end
end
