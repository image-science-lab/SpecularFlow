clear all
close all

addpath('./functions');
addpath('./AgrawalICCV05MatlabCode');

surfChoice = 42; %%%DEMO PARAMETER: Surface choice (20 = Quadric, 21 = Gaussian, 22 = Mixture of Gaussian, 23 = Buddha head)
EnvMap = 5; %%%Demo PARAMETER: Environment Value
RotationA = [ 0 5 10  ]; % 25 30 35 40 45 50]; %%Demo Parameter Rotataions
RotationShift = [ 5 ]; %%Demo Parameter Rotataions
interactive = 0; %DEMO PARAMETER Interactive Normals. Set this to "1" 
NormalsUnif = 0; %DEMO PARAMETER Uniform no of Normals. Set to 16 

pSize = 5; %%% DEMO PARAMETER: Patch size
throwBadFlow = 1; %%DEMO PARAMETER:

verb = 0; %%% DEMO PARAMETER: Verbose
saveResults = 0; %%DEMO PARAMETER


Rotation = [RotationA(:) RotationA(:)+RotationShift]';

if (1)
    createSyntheticImages(surfChoice, EnvMap, unique(Rotation(:)));
end

close all
[X_fl, Y_fl, f_fl, fx_fl, fy_fl, Zprime, fx_recons, fy_recons]= shapeFromSpecularFlow_SIFT(surfChoice, EnvMap, Rotation, pSize, throwBadFlow, NormalsUnif, interactive, verb, saveResults);


return
%%%%%%%%
if (0)
    sCount = 0;
    SNR = [];
    for surf=[20 21 22 41 42]
        %%Create Images
        if (0)
            createSyntheticImages(surf, EnvMap, unique(Rotation(:)));
        end

        if (0)
            for imC = EnvMap
                for rotC = 1:size(Rotation, 2)

                    img0 = imread( sprintf('../renderedImages/Surf%02d_Rot%02d_Env%02d.png', surf, Rotation(1, rotC), imC));
                    img1 = imread( sprintf('../renderedImages/Surf%02d_Rot%02d_Env%02d.png', surf, Rotation(2, rotC), imC));

                    for ch=1:3
                        %Compute SIFT Features
                        [F1, D1] = vl_sift(single(double(img0(:,:,ch))));
                        [F2, D2] = vl_sift(single(double(img1(:,:,ch))));

                        fIndx = (F1(3, :) > 10);
                        F1(:, fIndx) = []; D1(:, fIndx) = [];
                        fIndx = (F2(3, :) > 10);
                        F2(:, fIndx) = []; D2(:, fIndx) = [];


                        %Matching
                        [matc, sco] = vl_ubcmatch(D1, D2, 1.5);

                        FStruct{ch}.F1 = F1;            FStruct{ch}.F2 = F2;
                        FStruct{ch}.D1 = D1;            FStruct{ch}.D2 = D2;
                        FStruct{ch}.matc = matc;            FStruct{ch}.sco = sco;


                    end
                    save(sprintf('SIFT_Features/Surf%02d_Env%02d_RotStart%02d_RotStop%02d.mat', surf, imC, Rotation(1, rotC), Rotation(2, rotC)), 'FStruct');
                    clear FStruct
                end
            end
        end





        sCount = sCount + 1;
        for rotIndx = 1:size(Rotation, 2)
            Rot = Rotation(:, 1:rotIndx);
            [X_fl, Y_fl, f_fl, fx_fl, fy_fl, Zprime, fx_recons, fy_recons] = shapeFromSpecularFlow_SIFT_PreComputed(surf, EnvMap, Rot, pSize, throwBadFlow, 0, 0, 0, 0);

            SNR(sCount, rotIndx) = 10*log10(var(f_fl(:))/mean((f_fl(:)-Zprime(:)).^2));
            SNR;
            pause(1)
            close all

        end

    end



end