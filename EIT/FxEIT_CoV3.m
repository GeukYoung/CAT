function [CoVx, CoVy] = FxEIT_CoV3(Sigma,ROI,display_flag,maxscale)
%% prepare grid image
display_flag = 1;
maxscale = 4000;
% maxscale = 2550;
temp = size(Sigma);
temp(temp==1) = [];
Sigma_dim = length(temp);

Image = Sigma;

%% find FOV
if nargin < 2
    temp = ~isnan(Image);
else
    temp = ROI;
end
temp(~isnan(temp)) = 1;
FOV_h = nansum(temp,1);
FOV_h(FOV_h > 0) = 1;
nFOV_h = sum(FOV_h);
FOV_h1 = sum(~FOV_h(1:round(length(FOV_h)*0.5))) + 1;
FOV_h2 = FOV_h1 + nFOV_h - 1;

FOV_v = nansum(temp,2);
FOV_v(FOV_v > 0) = 1;
nFOV_v = sum(FOV_v);
FOV_v1 = sum(~FOV_v(1:round(length(FOV_v)*0.5))) + 1;
FOV_v2 = FOV_v1 + nFOV_v - 1;

FOV_Image = zeros(nFOV_v,nFOV_h);
FOV_Image = abs(Image(FOV_v1:FOV_v2,FOV_h1:FOV_h2));

%% CoV Calcuration
% CoVx
imagesc(FOV_Image);
total_sum = nansum(nansum(FOV_Image))
temp = 0;
for i = 1:size(FOV_Image,2);
    temp = temp + nansum(FOV_Image(:,i)) * (i/size(FOV_Image,2));
end
CoVx = temp/total_sum

% CoVy
total_sum = nansum(nansum(FOV_Image))
temp = 0;
for i = 1:size(FOV_Image,1);
    temp = temp + nansum(FOV_Image(i,:)) * (i/size(FOV_Image,1));
end
CoVy = temp/total_sum
disp(['CoVx : ',num2str(round(CoVx*1000)*0.1),'% , CoVy : ',num2str(round(CoVy*1000)*0.1),'%']);

raw_CoVx = CoVx * size(FOV_Image,2);
raw_CoVy = CoVy * size(FOV_Image,1);

%% Display result
% try
%     if display_flag == 1
%         subplot(2,3,[1 4]);
%         imagesc(Image); axis image off; colormap hot
%         subplot(2,3,[2 3 5 6]);
%         imagesc(-FOV_Image); axis image;
%         hold on;
%         plot([0 raw_CoVx-8],[raw_CoVy raw_CoVy],'r','LineWidth',1);
%         plot([raw_CoVx raw_CoVx],[0 raw_CoVy-8],'r','LineWidth',1);
%         plot(raw_CoVx,raw_CoVy,'r+','LineWidth',2); axis image off;
%         title(['CoVx : ',num2str(round(CoVx*1000)*0.1),'% , CoVy : ',num2str(round(CoVy*1000)*0.1),'%']);
%         set(gcf,'color','white');
%     end
% end

Cmap3 = FxImage_Cmap(3);
Cmap3 = Cmap3(end:-1:1,:); % reverse colormap
Cmap3(1,:) = [0 0 0]; % background white
try
    if display_flag == 1
        imagesc(Image, [-maxscale maxscale]); axis image off; colormap hot
        FxEIT_PlotBoundary(isnan(Image),'w',2);
        colormap(Cmap3);
        hold on;
        plot([FOV_h1 raw_CoVx+FOV_h1-15],[raw_CoVy+FOV_v1 raw_CoVy+FOV_v1],'r','LineWidth',2);
        plot([raw_CoVx+FOV_h1 raw_CoVx+FOV_h1],[FOV_v1 raw_CoVy+FOV_v1-15],'r','LineWidth',2);
        plot(raw_CoVx+FOV_h1,raw_CoVy+FOV_v1,'r+','LineWidth',2); axis image off;
%         title(['CoVx : ',num2str(round(CoVx*1000)*0.1),'% , CoVy : ',num2str(round(CoVy*1000)*0.1),'%']);
        
        % ROI Box
        plot([FOV_h1,FOV_h2],[FOV_v1,FOV_v1],'r','LineWidth',2);
        plot([FOV_h2,FOV_h2],[FOV_v1,FOV_v2],'r','LineWidth',2);
        plot([FOV_h1,FOV_h1],[FOV_v2,FOV_v1],'r','LineWidth',2);
        plot([FOV_h1,FOV_h2],[FOV_v2,FOV_v2],'r','LineWidth',2);
        FOV_Image = abs(Image(FOV_v1:FOV_v2,FOV_h1:FOV_h2));
        set(gcf,'color','white');
    end
end
%%
% % find ROI
% ROI_Image2 = abs(Image);
% ROI_Image2(ROI_Image2 < max(max(ROI_Image2))*th) = 0;
% ROI_Image2(isnan(ROI_Image2)) = 0;
% ROI_Image2(ROI_Image2 > 0) = 1;
% imagesc(ROI_Image2);
% 
% ROI_h = sum(ROI_Image2,1);
% ROI_h(ROI_h > 0) = 1;
% ROI_h1 = find(ROI_h,1,'first');
% ROI_h2 = find(ROI_h,1,'last');
% 
% ROI_v = sum(ROI_Image2,2)';
% ROI_v(ROI_v > 0) = 1;
% ROI_v1 = find(ROI_v,1,'first');
% ROI_v2 = find(ROI_v,1,'last');
% 
% ROI_Image = zeros(ROI_v2-ROI_v1+1,ROI_h2-ROI_h1+1);
% ROI_Image = Image(ROI_v1:ROI_v2,ROI_h1:ROI_h2);
% imagesc(ROI_Image); axis image off;




