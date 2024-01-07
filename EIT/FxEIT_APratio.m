function [AP_ratio] = FxEIT_APratio(Sigma,ROI,display_flag,maxscale)
%% prepare grid image
temp = size(Sigma);
temp(temp==1) = [];
Sigma_dim = length(temp);

Image = Sigma;
%% find FOV
temp = ~isnan(Image);
FOV_h = sum(temp,1);
FOV_h(FOV_h > 0) = 1;
nFOV_h = sum(FOV_h);
FOV_h1 = sum(~FOV_h(1:round(length(FOV_h)*0.5))) + 1;
FOV_h2 = FOV_h1 + nFOV_h - 1;

FOV_v = sum(temp,2);
FOV_v(FOV_v > 0) = 1;
nFOV_v = sum(FOV_v);
FOV_v1 = sum(~FOV_v(1:round(length(FOV_v)*0.5))) + 1;
FOV_v2 = FOV_v1 + nFOV_v - 1;

FOV_Image = zeros(nFOV_v,nFOV_h);
FOV_Image = Image(FOV_v1:FOV_v2,FOV_h1:FOV_h2);

FOV_half = round(nFOV_v * 0.5);
A_value = abs(nansum(nansum(FOV_Image(1:FOV_half,:))));
P_value = abs(nansum(nansum(FOV_Image(FOV_half+1:end,:))));
AP_ratio = A_value/P_value;
disp(['AP : ' ,num2str(AP_ratio)]);

%% Display result
try
    if display_flag == 1
%         subplot(221);
%         imagesc(Image); axis image off; colormap hot
%         subplot(223);
%         imagesc(FOV_Image); axis image off;
%         subplot(2,2,[2 4]);
        normA = A_value/(A_value + P_value);
        normP = P_value/(A_value + P_value);
        temp = [0 0;0 0.5;normP 0.5;normP 0;0 0.5;0 1;normA 1;normA 0.5];
%         patch('Vertices',temp,'faces',[1 2 3 4;5 6 7 8],'FaceColor',[0 0 0]); axis equal off; axis([0 1 0 1]);
%         hold on; plot([0 1],[0.5 0.5],'r','LineWidth',2);
%         plot([0 1],[0 0],'r','LineWidth',2);
%         plot([0 1],[1 1],'r','LineWidth',2);
%         plot([0 0],[0 1],'r','LineWidth',2);
%         plot([1 1],[0 1],'r','LineWidth',2);
% %         text(0.05, 0.25,'\Sigma P','fontsize',20,'Color','w');
% %         text(0.05, 0.75,'\Sigma A','fontsize',20,'Color','w');
%         text(0.80, 0.25,'P','fontsize',20,'Color','k');
%         text(0.80, 0.75,'A','fontsize',20,'Color','k');
%         set(gcf,'color','white');
%         title(['AP : ' ,num2str(AP_ratio)]);
        patch('Vertices',temp,'faces',[1 2 3 4;5 6 7 8],'FaceColor',[0.0779428571428572	0.503985714285714	0.838371428571429]); axis equal off; axis([0 1 0 1]);
        hold on; plot([0 1],[0.5 0.5],'w','LineWidth',2);
        plot([0 1],[0 0],'w','LineWidth',2);
        plot([0 1],[1 1],'w','LineWidth',2);
        plot([0 0],[0 1],'w','LineWidth',2);
        plot([1 1],[0 1],'w','LineWidth',2);
%         text(0.05, 0.25,'\Sigma P','fontsize',20,'Color','w');
%         text(0.05, 0.75,'\Sigma A','fontsize',20,'Color','w');
        text(0.80, 0.25,'P','fontsize',20,'Color','w');
        text(0.80, 0.75,'A','fontsize',20,'Color','w');
        set(gcf,'color','black');
        title(['AP : ' ,num2str(AP_ratio)]);
    end
end

%%
