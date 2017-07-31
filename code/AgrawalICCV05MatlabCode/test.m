
clear all;close all;clc;


ADD_NOISE = 0
ADD_OUTLIERS = 0

% Generate a synthetic surface im
synthetic_ramppeaks

im = im - min(im(:));

%find gradients
[ogx,ogy] = calculate_gradients(im,0,0);

% add noise
if(ADD_NOISE)
    
    % sigma = 5 percent of max gradient magnitude
    tt = sqrt(ogx.^2 + ogy.^2);
    sigma = 5*max(tt(:))/100
    clear tt
else
    sigma = 0
end


gx = ogx + sigma*randn(H,W);
gy = ogy + sigma*randn(H,W);


%add outliers
if(ADD_OUTLIERS)

    %find random locations to add outliers
    fac = 3
    outlier_x = rand(H,W) > 0.90;
    outlier_x(:,end) = 0;

    gx = gx + fac*outlier_x.*(2*(rand(H,W)>0.5)-1);

    outlier_y = rand(H,W) > 0.90;
    outlier_y(end,:) = 0;

    gy = gy + fac*outlier_y.*(2*(rand(H,W)>0.5)-1);
    outlier_x = double(outlier_x);
    outlier_y = double(outlier_y);
    disp(sprintf('Gx outliers = %d',sum(outlier_x(:))));
    disp(sprintf('Gy outliers = %d',sum(outlier_y(:))));
end










%Wflag = 0  use Least squares;
%Wflag = 1 if using abs curl as weights in ICCV algo
%Wflag = 2 if using grad mag as weights in ICCV algo
% else using abs(gx) abs(gy) as weights in ICCV algo
Wflag =2

% Boundary image: If given, then uses Dirichlet, else uses Neumann
bimage = []

% Threshold for deciding which gradients are in error
% If not specified use 10% of max curl value
TH_CURL = 1


Z = integrability_2d(gx,gy,Wflag,TH_CURL,bimage);


mydisplay(im);title('initial surface');
mydisplay(Z);title('reconstructed surface');









