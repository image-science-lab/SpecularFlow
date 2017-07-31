run DEMO.m




SYMBOLS used in various code.

xyprofile		Analytical depth map descriptor

f_depth  		ground truth depth map (high res)
f_grad_x f_grad_y	ground truth gradients (high res)
xGrid yGrid		X-Y grid for above
xRes yRes		X-Y resolution (scalars)

eigH eigL		Hessian eigenvalues high and low for surface
xFl yFl			X-Y grid for above

R Rt			Rotation values



f_fl			Depth map at flow's resolution
fx_fl fy_fl		Gradients at flow's resolution
X_fl Y_fl		X-Y grid for above


flow_x flow_y		Specular Flow displacements
mScore			Angular error in Flow in degrees

pSizeR pSizeC		Patch Size of Quadric approx
pShiftR pShiftC		When overlapping, this is the "shift"


OvMat1			Neigborhood constraints
OvMat2			Flow Constraints
OvMat3			Regularization Constraints (H = 0)
OvMat			Stacking all the above constraints
OvMat_Sq		full(OvMat'*OvMat)


xAns			Estimated Surface parameters
fx_recons fy_recons	Reconstructed gradients
Zprime			Reconstructed surface

X_sta Y_sta		SIFT Point First image
X_stp Y_stp		SIFT points second images
dis			Accuracy of SIFT Matches

