Codebase for 

Sankaranarayanan, A. C.; Veeraraghavan, A.; Tuzel, O.; and Agrawal, A. Specular surface reconstruction from sparse reflection correspondences. In IEEE Intl. Conf. Computer Vision and Pattern Recognition (CVPR), 2010

@inproceedings{sankaranarayanan2010specular,
  title={Specular surface reconstruction from sparse reflection correspondences},
  author={Sankaranarayanan, Aswin C. and Veeraraghavan, Ashok and Tuzel, Oncel and Agrawal, Amit},
  booktitle={IEEE Intl. Conf. Computer Vision and Pattern Recognition (CVPR)},
  year={2010},
  url_Paper = {files/paper/2010/SfRC_CVPR2010.pdf},
  pubtype = {Conference},
}


YOU NEED CVX AND VL_FEAT to run the code in this package. get them first!!

then run DEMO.m in code folder

images folder has some environment textures.

renderedimages are where the renderings of the mirror will be put.


Pre-use stuff:
1) Get CVX Configured. Do this by going to cvx folder and running cvx_setup inside MATLAB. it maybe worthwhile adding DEMO/cvx directory tree to ur MATLAB Path.

2) SIFT Toolbox: VLFEAT: To get it setup, you need to go to vlfeat/toolbox and run vl_setup.m . To test, u can run vl_demo.m [NOTE that this toolbox demo does not like blank spaces in the base directory structure. The SIFT part itself might work, but the demo will crash if your base directory has blank spaces]. Add the directory tree to your path.


If you dont want to save these paths, just run cvx_setup and vl_setup each time.
------------