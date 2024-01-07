function [Quad_UL, Quad_UR, Quad_LL, Quad_LR] = FxImage_Quadsplit(Image,op)
    if nargin > 1
%         [half_vertical,half_horizontal] = size(Image);
%         half_horizontal = round(half_horizontal*0.5);
%         half_vertical = round(half_vertical*0.5);
        half_horizontal = op(1);
        half_vertical = op(2);
    else
        figure; imagesc(Image); set(gca,'Ydir','Normal');
        background_img = Image;
        background_img(background_img==0) = 1;
        FxEIT_PlotBoundary(background_img,'w',3);
        temp = round(ginput(1));
        half_horizontal = temp(1);
        half_vertical = temp(2);
        close;
    end
%     Quad_LL = Image(1:half_vertical,1:half_horizontal);
%     Quad_LR = Image(1:half_vertical,half_horizontal+1:end);
%     Quad_UL = Image(half_vertical+1:end,1:half_horizontal);
%     Quad_UR = Image(half_vertical+1:end,half_horizontal+1:end);

%     figure;
%     subplot(2,4,[1 2 5 6]); imagesc(Image); axis image off;
%     subplot(2,4,3); imagesc(Quad_LL); axis image off;
%     subplot(2,4,4); imagesc(Quad_LR); axis image off;
%     subplot(2,4,7); imagesc(Quad_UL); axis image off;
%     subplot(2,4,8); imagesc(Quad_UR); axis image off;
    
    Quad_UL = Image;
    Quad_UL(1:half_vertical,:) = 0;
    Quad_UL(:,half_horizontal+1:end) = 0;
    
    Quad_UR = Image;
    Quad_UR(1:half_vertical,:) = 0;
    Quad_UR(:,1:half_horizontal) = 0;
    
    Quad_LL = Image;
    Quad_LL(half_vertical+1:end,:) = 0;
    Quad_LL(:,half_horizontal+1:end) = 0;
    
    Quad_LR = Image;
    Quad_LR(half_vertical+1:end,:) = 0;
    Quad_LR(:,1:half_horizontal) = 0;
    
    figure;
    subplot(2,4,[1 2 5 6]); imagesc(Image); axis image off; set(gca,'Ydir','Normal'); hold on;
    plot(half_horizontal,half_vertical,'ro'); FxEIT_PlotBoundary(background_img,'w',3);
    subplot(2,4,3); imagesc(Quad_UL); axis image off; set(gca,'Ydir','Normal'); title('UL(1)');set(gca,'Ydir','Normal');
    subplot(2,4,4); imagesc(Quad_UR); axis image off; set(gca,'Ydir','Normal'); title('UR(2)');set(gca,'Ydir','Normal');
    subplot(2,4,7); imagesc(Quad_LL); axis image off; set(gca,'Ydir','Normal'); title('LL(3)');set(gca,'Ydir','Normal');
    subplot(2,4,8); imagesc(Quad_LR); axis image off; set(gca,'Ydir','Normal'); title('LR(4)'); set(gca,'Ydir','Normal');

end