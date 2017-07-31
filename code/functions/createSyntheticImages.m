function  createSyntheticImages(surfChoice, imChoice, rotAngle)

disp('Synthesizing image ....');

%%%%%%%%%%%Analytical surface descriptions
%%%SURFACE CREATION
[xyprofile, xGrid, yGrid, cutoff, analytical, f_depth ] = surfaceCreator(surfChoice);



for rtC = 1:length(rotAngle)
    rt = rotAngle(rtC);
    omega = [0 0.01; -0.01 0]*rt;
    R = expm(omega); Rt = R';

    if (analytical)
        x = Rt(1,1)*xGrid + Rt(1,2)*yGrid;
        y = Rt(2,1)*xGrid + Rt(2,2)*yGrid;
        frot_depth = eval(xyprofile);
        if ~isempty(cutoff)
            frot_depth = (frot_depth >= cutoff).*(frot_depth - cutoff);
        end
    else
        x = Rt(1,1)*xGrid + Rt(1,2)*yGrid;
        y = Rt(2,1)*xGrid + Rt(2,2)*yGrid;
        frot_depth = interp2(xGrid, yGrid, f_depth, x, y, 'bicubic', 0);
    end


    for imCC = 1:length(imChoice)
        imC = imChoice(imCC);
        envMap = imread(sprintf('../images/%08d.png', imC)); envMap = double(envMap)/255;

        img1 = renderSurface_equiArea(frot_depth,  xGrid, yGrid, envMap);
        imshow([img1 ]);
        drawnow

        imwrite(img1, sprintf('../renderedImages/Surf%02d_Rot%02d_Env%02d.png', surfChoice, rt, imC));

    end
end