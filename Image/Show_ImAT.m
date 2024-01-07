function [Image_AT] = Show_ImAT(Image,DataSet,cmap,cmax)
% input parameters
% 1. Image - 1xP^2: reshape by sqrt(length)
%          - PxP
% 2. DataSet - using boundary image & colormap(when w/o cmap option) & RM(when Image reshape)
% 3. cmap - 1: oneside(B/P)
%         - 2: doubleside(B/R)
%         - Cx3: input colormap(suggest, length C over 256) 
% 4. cmax - none: scale fit
%         - abcmax: [-abcmax abcmax]
%         - [cmin cmax]: [cmin cmax]
%
% output parameter
% 1. Image - PxP(colormap + boundary shape)

if nargin < 4
    cmax = max(max(abs(Image)))*1.1;
end

if nargin < 3
    cmap = 1;
end

if min(size(cmap)) == 1
    switch cmap
        case 1
            cmap = DataSet.EIT.Image.Cmap_oneside;
        case 2
            cmap = DataSet.EIT.Image.Cmap_doubleside;
    end
elseif size(cmap,2) ~= 3
    cmap = cmap';
end

if size(Image,1) == 1 || size(Image,2) == 1
    if length(Image) == size(DataSet.EIT.Image.RM,2)
        Image = (DataSet.EIT.Image.RM * Image);
    end
    Image = reshape(Image, sqrt(size(DataSet.EIT.Image.RM,1)), sqrt(size(DataSet.EIT.Image.RM,1)));
    % AT <-> Matlab since x,y coordination diff
    Image = Image';
end

F = griddedInterpolant({1:64 1:64},Image);
F_intp = F({linspace(1,64,2048) linspace(1,64,2048)});
if length(cmax) == 1
    F_intp = uint8(((F_intp)/2/cmax+0.5)*255);
    Image_AT = ind2rgb(F_intp, cmap);
elseif length(cmax) == 2
    if cmax(1) > cmax(2)
        cmax = [cmax(2) cmax(1)];
    end
%     idx_none = F_intp<cmax(1);
    F_intp = uint8(((F_intp-cmax(1))/diff(cmax))*256);
%     F_intp(idx_none) = 0;
    Image_AT = ind2rgb(F_intp, cmap);
%     Image_AT(:,:,1) = Image_AT(:,:,1).*~idx_none;
%     Image_AT(:,:,2) = Image_AT(:,:,2).*~idx_none;
%     Image_AT(:,:,3) = Image_AT(:,:,3).*~idx_none;
end

Image_AT(:,:,1) = Image_AT(:,:,1).*(DataSet.EIT.Image.mask_alpha==0);
Image_AT(:,:,2) = Image_AT(:,:,2).*(DataSet.EIT.Image.mask_alpha==0);
Image_AT(:,:,3) = Image_AT(:,:,3).*(DataSet.EIT.Image.mask_alpha==0);
Image_AT = Image_AT + double(DataSet.EIT.Image.mask_cdata)/256;

% figure; imagesc(Image_AT); axis image off;
% figure; imagesc(F_intp); colormap('copper')
% figure; plot(F_intp(:)); 
% figure; plot(Image_AT(:)); 
% figure; imagesc(F_intp); colormap copper;
end

