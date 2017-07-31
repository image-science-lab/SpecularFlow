

Matlab code

For research purposes only. If you are using this code for a journal or conference publication, please include the following citations

A. Agrawal, R. Chellappa and R. Raskar, "An algebraic approach to surface reconstruction from gradient fields". International conference on computer vision (ICCV) 2005

A. Agrawal, R. Raskar, and R. Chellappa, "What is the range of surface reconstructions from a gradient field?" European conference on Computer Vision (ECCV) 2006


For feedback and suggestions, please email  aagrawal at umd dot edu

Copyright (c) Amit Agrawal 2004-2006 
http://www.cfar.umd.edu/~aagrawal


At matlab prompt, run

-> test

This will generate a synthetic surface, find its gradients, add noise and outliers to it and reconstruct the surface.


Notes:

-If you just want get a least square solution, use Wflag = 0. This finds a solution by solving the Poisson equation. My implementation is exact analytically and does not require any magic parameters. 

-If the gradient field has no noise and no outliers, this implementation gives an exact solution up to an additive constant. Note that techniques based on multigrid always have residual errors even in this case.




-For large grid size, the algorithm may take time to run as it finds the Minimum spanning tree (MST) of the underlying graph. The MST code is written in C, is optimized and included in Matlab via dll.


