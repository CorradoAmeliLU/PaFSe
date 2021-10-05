function F = getFeatFast(I,SE,GPU)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Praqa is an algorithm for quantifying protein relative abundance from
% 3D fluorescent acquisitions. The algorithm reproduces the results from
% the article:
% PRAQA: Protein Relative Abundance Quantification Algorithm for 3D Fluorescent Images
% Corrado Ameli, Sonja Fixemer, David S. Bouvier, Alexander Skupin
%
%
% INPUT:
% - I : The 3D image in double format
% - SE : List of SEs
% - CleanProc : Enables or disables cleaning of isolated positive signal 
%               pixels.
% - GPU: Uses GPU for processing image.
%
%
% OUTPUT:
% - F : [Numer of Pixels BY Number of SEs] matrix containing average
%       intensity values for each pixel for a particular SE.
%
%
% All Rights Reserved
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

se = SE.Neighborhood;
n_ones = sum(se==1,'all');
se = double(se)./n_ones;

if GPU
    I = gpuArray(I);
    se = gpuArray(se);
end

F = conv2(I,se,'valid');

[Fx,Fy] = size(F);
[Ix,Iy] = size(I);

PADx = (Ix - Fx)/2;
PADy = (Iy - Fy)/2;

if GPU
    F_pad = nan(size(I),'gpuArray');
else
    F_pad = nan(size(I));
end

F_pad(PADx+1 : Ix-PADx , PADy+1 : Iy-PADy) = F;
F = reshape(F_pad,[Ix*Iy 1]);

if GPU
    F = gather(F);
end

end