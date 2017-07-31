function [flow_x, flow_y, matchScore] = computeSpecularFlow(xGrid, yGrid, f_depth, frot_depth, R,  X_fl, Y_fl)



xRes = xGrid(1,2)-xGrid(1,1);
yRes = yGrid(2,1)-yGrid(1,1);


[f_grad_x, f_grad_y] = gradient(f_depth, xRes, yRes);
[fnew_grad_x, fnew_grad_y] = gradient(frot_depth, xRes, yRes);

%X_fl = xGrid; fx_fl = f_grad_x;
%Y_fl = yGrid; fy_fl = f_grad_y;

fx_fl = interp2(xGrid, yGrid, f_grad_x, X_fl, Y_fl);
fy_fl = interp2(xGrid, yGrid, f_grad_y, X_fl, Y_fl);



flow_x = zeros(size(X_fl));
flow_y = zeros(size(X_fl));


for row=1:size(X_fl, 1)

    disp(size(X_fl,1)-row)
    for col=1:size(X_fl, 2)


        idx = find(abs(xGrid(1,:) - X_fl(row, col)) <  xRes/2);
        idy = find(abs(yGrid(:,1) - Y_fl(row, col)) <  yRes/2);
        
        angMatch  = 1e10;
        if (f_depth(idy, idx) == 0)
            local_flow = 0/0*ones(2,1);
        else

            siz = 5; shft = 1;
            local_fx = f_grad_x(max(1,idy-siz):shft:min(size(f_grad_x, 1),idy+siz), max(1,idx-siz):shft:min(size(f_grad_x, 2), idx+siz));
            local_fy = f_grad_y(max(1,idy-siz):shft:min(size(f_grad_x, 1),idy+siz), max(1,idx-siz):shft:min(size(f_grad_x, 2), idx+siz));
            local_x = xGrid(max(1,idy-siz):shft:min(size(xGrid, 1),idy+siz), max(1,idx-siz):shft:min(size(xGrid, 2), idx+siz));
            local_y = yGrid(max(1,idy-siz):shft:min(size(xGrid, 1),idy+siz), max(1,idx-siz):shft:min(size(xGrid, 2), idx+siz));

            Amat = [ 1+0*local_x(:) 0*local_x(:) local_x(:) local_y(:)  0*local_x(:);
                0*local_x(:) 1+0*local_x(:) 0*local_x(:) local_x(:) local_y(:);];
            Bmat = [ local_fx(:); local_fy(:)];

            param = Amat\Bmat;

            Jx0 = [param(1); param(2)];
            Hx0 = [ param(3) param(4); param(4) param(5)];

            if sum(abs(Hx0)) == 0
                local_flow = 0/0 * ones(2,1);
            else
                local_flow = -(eye(2)-R*pinv(Hx0)*R'*Hx0)*[X_fl(row,col); Y_fl(row,col)]- R*pinv(Hx0)*(eye(2)-R')*Jx0;
            end

            siz0 = 30;

            idy1 = (idy+round(local_flow(2)/xRes));
            idx1 = (idx+round(local_flow(1)/xRes));

            if (idy1 > 1) & (idy1 < size(xGrid, 1)) &(idx1 > 1) & (idx1 < size(xGrid,2))

                fx_tmp1 = fnew_grad_x(max(1, idy1-siz0):min(size(xGrid, 1), idy1+siz0), max(1, idx1-siz0):min(size(xGrid, 2), idx1+siz0));
                fy_tmp1 = fnew_grad_y(max(1, idy1-siz0):min(size(xGrid, 1), idy1+siz0), max(1, idx1-siz0):min(size(xGrid, 2), idx1+siz0));
                x_tmp1 = xGrid(max(1, idy1-siz0):min(size(xGrid, 1), idy1+siz0), max(1, idx1-siz0):min(size(xGrid, 2), idx1+siz0));
                y_tmp1 = yGrid(max(1, idy1-siz0):min(size(xGrid, 1), idy1+siz0), max(1, idx1-siz0):min(size(xGrid, 2), idx1+siz0));

                res = 1/10;
                [x_tmp, y_tmp] = meshgrid(min(x_tmp1(:)):xRes*res:max(x_tmp1(:)), min(y_tmp1(:)):yRes*res:max(y_tmp1(:)));
                fx_tmp = interp2(x_tmp1, y_tmp1, fx_tmp1, x_tmp, y_tmp, 'spline');
                fy_tmp = interp2(x_tmp1, y_tmp1, fy_tmp1, x_tmp, y_tmp, 'spline');

                dMap = (fx_fl(row, col) - fx_tmp).^2 + (fy_fl(row, col) - fy_tmp).^2;

                [maxVal] = min(dMap(:));
                indx = find(dMap == maxVal);
                [irow, icol] = ind2sub(size(dMap), indx);

                zz = dMap < maxVal*10;
                [zz1, num] = bwlabel(zz);
                if (num > 1)
                    for qq=1:num
                        tmpD = dMap.*(zz1 == qq) + 2000000*(zz1 ~= qq);
                        [kk1, kk2] = min(tmpD(:));
                        [irow1, icol1] = ind2sub(size(dMap), kk2(1));
                        if (abs(irow1-siz0/res-1)^2+abs(icol1-siz0/res-1)^2) < (abs(irow-siz0/res-1)^2+abs(icol-siz0/res-1)^2)
                            irow = irow1;
                            icol = icol1;
                        end
                    end
                end

                if (length(irow)  == 1)
                    %%%%MinMatchStuff
                    angMatch = (fx_fl(row, col)*fx_tmp(irow, icol)+fy_fl(row,col)*fy_tmp(irow, icol)+1)/(norm([-fx_fl(row,col); -fy_fl(row,col);1])*norm([-fx_tmp(irow, icol); -fy_tmp(irow, icol); 1]));
                    angMatch = acos(angMatch)*180/pi;

                    local_flow = [x_tmp(irow, icol)-X_fl(row, col); y_tmp(irow, icol)-Y_fl(row, col)]; %local_flow + [icol-siz0-1;irow-siz0-1]*xRes;
                else
                    local_flow = 0/0*ones(2,1);
                end

                if (1)
                    subplot 121
                    imagesc(flow_x~= 0);
                    subplot 122
                    hold off
                    imagesc(log10(dMap));
                    hold on
                    plot(siz0/res+1, siz0/res+1, 'r*');
                    plot(icol, irow, 'ro');
                    drawnow
                end

            else
                local_flow = 0/0*ones(2,1);
            end
        end


        flow_x(row, col) = local_flow(1);
        flow_y(row, col) = local_flow(2);
        matchScore(row, col) = angMatch;
    end
end