%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PaFSe is an algorithm for segmenting 3D fluorescent images
% The algorithm reproduces the experimental
% results from the article:
%
% PaFSe: a Parameter Free Segmentation Approach for 
%        3D Fluorescent Images (Ameli et al)
%  
%
% INPUT:
% - I :             The 3D image in double format
% - parallelize:    If set to true, multiple slices will be processed
%                   at the same time.
% - disp_progress : Display progress in the console.
%
%
% OUTPUT:
% - I_out : Segmented 3D image.
% - vol : Count of all the positive pixels of I_out.
% - est_thetas: Theta parameters estimated for each slice.
%
%
% All Rights Reserved
%
%
% REQUIRES IMAGE PROCESSING TOOLBOX and PARALLEL COMPUTING TOOLBOX
% ONLY IF parallelize = TRUE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [I_out,vol,est_thetas] = PaFSe(I,parallelize,disp_progress)
    
    addpath(genpath(('./Library/')));

    [X,Y,Z] = size(I);
    I_out = zeros(X,Y,Z);

    est_thetas = nan(1,Z);

    if parallelize

        parfor z = 1 : Z
            
            fitEta = load('./Library/fitEta.mat');
            fitEta = fitEta.fitEta;
            
            if disp_progress
                disp(z)
            end

            i_slice = I(:,:,z);

            estimated_theta = round(feval(fitEta,estimate_noise(i_slice))*10,0);

            [~,I_out(:,:,z)] = Praqa(i_slice,estimated_theta/10,0,0);
            est_thetas(z) = estimated_theta/10;
        end
    else
        for z = 1 : Z
            
            load('./Library/fitEta.mat')
            
            if disp_progress
                disp(z)
            end

            i_slice = I(:,:,z);

            estimated_theta = round(feval(fitEta,estimate_noise(i_slice))*10,0);

            [~,I_out(:,:,z)] = Praqa(i_slice,estimated_theta/10,0,0);
            est_thetas(z) = estimated_theta/10;
        end

    end

    vol = sum(I_out,'all');

end

