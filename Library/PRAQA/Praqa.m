function [Volume,I_final] = Praqa(I,ThresholdFactor,CleanProc,GPU)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Praqa is an algorithm for quantifying protein relative abundance from
% 3D fluorescent acquisitions. The algorithm reproduces the experimental
% results from the article:
%
% PRAQA: Protein Relative Abundance Quantification Algorithm for 3D 
%        Fluorescent Images (Ameli et al)
%  
%
% INPUT:
% - I : The 3D image in double format
% - ThresholdFactor : Threshold for detecting positive signal (see Theta 
%                     in the manuscript).
% - CleanProc : Enables or disables cleaning of isolated positive signal 
%               pixels.
% - GPU: Uses GPU for processing image.
% - SE (optional) : to change the SEs, edit SE.mat file and add the 
%                   desired SEs.
%
%
% OUTPUT:
% - Volume : Sum of the positive pixels throughout the whole 3D image.
% - I_final : Binarized image of I.
%
%
% All Rights Reserved
%
%
% REQUIRES IMAGE PROCESSING TOOLBOX and PARALLEL COMPUTING TOOLBOX
% ONLY IF CLEANPROC = TRUE or GPU = TRUE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Ix,Iy,Iz] = size(I);

%RANGE (0,1) NORMALIZATION

MAX = max(I,[],'all');
I_norm = I./MAX;

%LOAD STRUCTURING ELEMENTS

SeList = load('SE.mat');
SeList = SeList.SE;
n_SE = length(SeList);

%CREATE FEATURE SET AND DETECT MARKER (SLICE WISE)

I_marked = zeros(Ix,Iy,Iz);

for zz = 1 : Iz
        
    feats = zeros(Ix*Iy,n_SE);
    I_filt = I_norm(:,:,zz);
    
    for ii = 1 : n_SE
        feats(:,ii) = getFeatFast(I_filt,SeList{ii},GPU);
    end
    
    isMarker = false(size(feats));
    for ii = 1 : n_SE
        [~,isMarker(:,ii)] = rmoutliers(feats(:,ii),'ThresholdFactor',ThresholdFactor);
    end
    
    isMarker = sum(isMarker,2)>0;
    
    I_fix = I_filt;
    I_fix(~reshape(isMarker,[Ix Iy])) = 0;
    
    I_marked(:,:,zz) = I_fix;
    
end

%CLEANING IMAGE FROM ISOLATED ARTIFACTS

if CleanProc
    
    I_marked(I_marked>0) = 1;
    I_cleared = logical(I_marked);
    
    for zz = 1 : Iz
        I_cleared(:,:,zz) = imopen(I_marked(:,:,zz),strel('square',2));
    end
    
else
    I_cleared = I_marked;
end

I_final = I_cleared > 0;

Volume = sum(I_final,'all');

end

