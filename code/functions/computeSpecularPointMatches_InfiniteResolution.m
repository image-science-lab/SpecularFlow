function [flow_x, flow_y, matchScore, flow_xi, flow_yi, iScore] = computeSpecularPointMatches(xGrid, yGrid, f_depth, frot_depth, xyprofile, R,  X_fl, Y_fl)

Rt = R';

xRes = xGrid(1,2)-xGrid(1,1);
yRes = yGrid(2,1)-yGrid(1,1);


[f_grad_x, f_grad_y] = gradient(f_depth, xRes, yRes);
[fnew_grad_x, fnew_grad_y] = gradient(frot_depth, xRes, yRes);

%X_fl = xGrid; fx_fl = f_grad_x;
%Y_fl = yGrid; fy_fl = f_grad_y;

fx_fl = interp2(xGrid, yGrid, f_grad_x, X_fl, Y_fl, 'bilinear');
fy_fl = interp2(xGrid, yGrid, f_grad_y, X_fl, Y_fl, 'bilinear');

clear x
clear y

syms x y m n
x = Rt(1,1)*m + Rt(1,2)*n;
y = Rt(2,1)*m + Rt(2,2)*n;

fmap = eval(xyprofile);
fxmap = diff(fmap, m);
fymap = diff(fmap, n);



flow_x = zeros(size(X_fl));
flow_y = zeros(size(X_fl));
flow_xi = zeros(size(X_fl));
flow_yi = zeros(size(X_fl));
matchScore = zeros(size(X_fl));
iScore = zeros(size(X_fl));


%%%normalized flow matching
normfl_1 = -fx_fl./sqrt(1+fx_fl.^2+fy_fl.^2);
normfl_2 = -fy_fl./sqrt(1+fx_fl.^2+fy_fl.^2);
normfl_3 = 1./sqrt(1+fx_fl.^2+fy_fl.^2);

normrot_1 = -fnew_grad_x./sqrt(1+fnew_grad_x.^2 + fnew_grad_y.^2);
normrot_2 = -fnew_grad_y./sqrt(1+fnew_grad_x.^2 + fnew_grad_y.^2);
normrot_3 = 1./sqrt(1+fnew_grad_x.^2 + fnew_grad_y.^2);


for row=1:size(X_fl, 1)

    disp(size(X_fl,1)-row)
    for col=1:size(X_fl, 2)


        idx = find(abs(xGrid(1,:) - X_fl(row, col)) <  xRes/2);
        idy = find(abs(yGrid(:,1) - Y_fl(row, col)) <  yRes/2);

        angMatch  = 0/0;        angMatch_i  = 0/0;
        local_flow = 0/0*ones(2,1);
        local_flow_i = 0/0*ones(2,1);

        %         if (f_depth(idy, idx) == 0)
        %             local_flow = 0/0*ones(2,1);
        %             local_flow_i = 0/0*ones(2,1);
        %         else

        angCost = (normfl_1(row, col)*normrot_1 + normfl_2(row, col)*normrot_2 + normfl_3(row, col)*normrot_3);

        [mVal, mIndx] = max(angCost(:));

        mVal = acos(mVal)*180/pi;

        [iy, ix] = ind2sub(size(normrot_1), mIndx);
        local_flow(1) = xGrid(iy, ix) - X_fl(row, col);
        local_flow(2) = yGrid(iy, ix) - Y_fl(row, col);
        angMatch = mVal;

        opt = optimset('TolFun', 1e-9, 'TolFun', 1e-9, 'TolX', 1e-9);
%        [X0, F0] = fsolve( @(a)   subs(subs([fxmap - fx_fl(row, col); fymap - fy_fl(row,col)], m, a(1)), n, a(2)), [ xGrid(iy, ix); yGrid(iy, ix)], opt);
        [X0, F0] = lsqnonlin( @(a)   subs(subs([fxmap - fx_fl(row, col); fymap - fy_fl(row,col)], m, a(1)), n, a(2)), [ xGrid(iy, ix); yGrid(iy, ix)], [], [], opt);



        local_flow_i(1) = X0(1) - X_fl(row,col);
        local_flow_i(2) = X0(2) - Y_fl(row,col);

        zx = [-subs(subs([fxmap; fymap], m, X0(1)), n, X0(2)); 1];
        zx =zx/norm(zx);
        [normtmp1] = zx(1);
        [normtmp2] = zx(2);
        [normtmp3] = zx(3);

        angMatch_i = acos((normfl_1(row, col)*normtmp1 + normfl_2(row, col)*normtmp2 + normfl_3(row, col)*normtmp3))*180/pi;




        %%%run small loop



        %             dMap = (fx_fl(row, col) - fnew_grad_x).^2 + (fy_fl(row, col) - fnew_grad_y).^2;
        %
        %             [maxVal] = min(dMap(:));
        %             indx = find(dMap == maxVal);
        %             [irow, icol] = ind2sub(size(dMap), indx);
        %
        %             irow = irow(1);
        %             icol = icol(1);
        %
        %             angMatch = (fx_fl(row, col)*fnew_grad_x(irow, icol)+fy_fl(row,col)*fnew_grad_y(irow, icol)+1)/(norm([-fx_fl(row,col); -fy_fl(row,col);1])*norm([-fnew_grad_x(irow, icol); -fnew_grad_y(irow, icol); 1]));
        %             angMatch = acos(angMatch)*180/pi;
        %
        %             local_flow = [xGrid(irow, icol)-X_fl(row, col); yGrid(irow, icol)-Y_fl(row, col)]; %local_flow + [icol-siz0-1;irow-siz0-1]*xRes;


        %         end


        flow_x(row, col) = local_flow(1);
        flow_y(row, col) = local_flow(2);
        matchScore(row, col) = angMatch;
        flow_xi(row, col) = local_flow_i(1);
        flow_yi(row, col) = local_flow_i(2);
        iScore(row, col) = angMatch_i;
    end
end