function [] = FxEIT_ComplianceCOD(PEEP,Image,ROI)
    % generate mask info (background & ROI) 
    Image(Image<0) = 0;
    background_mask = isnan(Image(:,:,1));
    background_mask_tab = reshape(background_mask,size(Image,1)*size(Image,2),1);
    ROI(isnan(ROI)) = 0;
    ROI_mask = logical(reshape(ROI,size(Image,1)*size(Image,2),1));
    
    % data tabulation
    for cntImage = 1:size(Image,3)
        Image_tab(:,cntImage) = reshape(Image(:,:,cntImage).*ROI,size(Image,1)*size(Image,2),1);
    end
    
    % find maximum compliance & pressure
    for cntPixel = 1:size(Image_tab,1)
        [Cmax_tab(cntPixel,1) Pstar_tab(cntPixel,1)] = max(Image_tab(cntPixel,:));
    end
    Pstar_tab(isnan(Pstar_tab)) = NaN;
    Pstar_tab = PEEP(Pstar_tab)';
    
    % calcurate collapse & overdistantion area
    for cntImage = 1:size(Image,3)
        for cntPixel = 1:size(Image_tab,1)
            if PEEP(cntImage) < Pstar_tab(cntPixel)
                Collapse_tab(cntPixel,cntImage) = ...
                    (Cmax_tab(cntPixel) - Image_tab(cntPixel,cntImage))/Cmax_tab(cntPixel);
                Overdist_tab(cntPixel,cntImage) = 0;
            elseif PEEP(cntImage) > Pstar_tab(cntPixel)
                Collapse_tab(cntPixel,cntImage) = 0;
                Overdist_tab(cntPixel,cntImage) = ...
                    (Cmax_tab(cntPixel) - Image_tab(cntPixel,cntImage))/Cmax_tab(cntPixel);
            else
                Collapse_tab(cntPixel,cntImage) = 0;
                Overdist_tab(cntPixel,cntImage) = 0;
            end
        end
        Collapse_tab(background_mask_tab,cntImage) = NaN;
        Overdist_tab(background_mask_tab,cntImage) = NaN;
        Collapse_tab(isinf(Collapse_tab)) = 0;
        Overdist_tab(isinf(Overdist_tab)) = 0;
        
        CollapseSum(cntImage) = round(nansum(Collapse_tab(:,cntImage).*Cmax_tab)/nansum(Cmax_tab)*100*10)/10;
        OverdistSum(cntImage) = round(nansum(Overdist_tab(:,cntImage).*Cmax_tab)/nansum(Cmax_tab)*100*10)/10; 
    end
    
    COD_tab = Overdist_tab - Collapse_tab;
    clim_scale = nanmax(nanmax(abs(COD_tab)))*1.01;
%     clim_scale = 1;
    
    % Set background and w/o ROI color using inf & NaN
    COD_tab(~ROI_mask,:) = inf;
    COD_tab(background_mask_tab,:) = inf;
    
    Collapse_tab(~ROI_mask,:) = inf;
    Collapse_tab(background_mask_tab,:) = inf;
    
    Overdist_tab(~ROI_mask,:) = inf;
    Overdist_tab(background_mask_tab,:) = inf;
        
    Pstar_tab(~ROI_mask) = inf;
    Pstar_tab(background_mask_tab) = inf;
    
    Cmax_tab(~ROI_mask) = inf;
    Cmax_tab(background_mask_tab) = inf;
    
    % convert tabulation to Image format
    for i = 1:size(Image,3)
        COD(:,:,i) = reshape(COD_tab(:,i),256,256);
        Collapse(:,:,i) = reshape(Collapse_tab(:,i),256,256);
        Overdist(:,:,i) = reshape(Overdist_tab(:,i),256,256);
    end
    Pstar(:,:) = reshape(Pstar_tab,256,256);
    Cmax(:,:) = reshape(Cmax_tab,256,256);
    
%     figure;
%     subplot(121);
%     imagesc(reshape(Cmax,256,256)); caxis([-50 50]); set(gca,'Ydir','normal');
%     subplot(122);
%     imagesc(reshape(Pstar,256,256)); set(gca,'Ydir','normal'); caxis([-30 30]);

    % 
%     Cmap3 = FxImage_Cmap(3);
%     Cmap3 = Cmap3(end:-1:1,:);
%     Cmap3(1,:) = [1 1 1];
%     Cmap4 = FxImage_Cmap(4);
%     for i = 1:size(Image,3)
%         subplot(3,size(Image,3),i);
%         imagesc(Collapse(:,:,i)); colormap(Cmap3); caxis([-1 1]); 
%         set(gca,'Ydir','normal'); axis image off;
%         
%         subplot(3,size(Image,3),i + size(Image,3));
%         imagesc(Overdist(:,:,i)); colormap(Cmap3); caxis([-1 1]); 
%         set(gca,'Ydir','normal'); axis image off;
%         
%         subplot(3,size(Image,3),i + 2*size(Image,3));
%         imagesc(COD(:,:,i)); colormap(Cmap3); caxis([-1 1]); 
%         set(gca,'Ydir','normal'); axis image off;
%     end
%     imagesc(Cmax); colormap(Cmap3); caxis([-120 120]); set(gca,'Ydir','normal'); axis image off;
%     imagesc(Pstar); colormap(Cmap3); caxis([-15 15]); set(gca,'Ydir','normal'); axis image off;
    
    %% Display PEEP COD Images
    Cmap3 = FxImage_Cmap(3); % import double side colormap
    Cmap3 = Cmap3(end:-1:1,:); % reverse colormap
    Cmap3(1,:) = [1 1 1]; % background white
%     Cmap3(end,:) = [];
    
    figure; set(gcf,'units','normalized','outerposition',[0.3 0.2 0.58 0.63]);
    whitebg('k'); set(gcf,'color','k');
    for i = 1:size(Image,3)
        subplottight(4,size(Image,3)+1,i);
        imagesc(COD(:,:,i)); caxis([-clim_scale clim_scale]); colormap(Cmap3); freezeColors;
%         image(COD(:,:,i));
%         set(gca,'Ydir','reverse'); 
        set(gca,'Ydir','normal'); 
        axis image off;
%         FxEIT_PlotBoundary(background_mask,'w',2);
        FxEIT_PlotBoundary(ROI,'w',1);
        title(['PEEP : ',num2str(PEEP(i))],'Color','w');
        set(gca, 'XTickLabel', [],'XTick',[])
    end
    subplottight(4,size(Image,3)+1,size(Image,3)+1);
    caxis([-clim_scale clim_scale]); colormap(Cmap3); colorbar('Ticks',[-clim_scale,clim_scale],'TickLabels',{'(Collapse)','(Overdist)'},'Location','west'); 
    axis off;  

    % Display Maximum Compliance Image
    subplottight(4,size(Image,3)+1,size(Image,3)+2);
    % colormap(Cmap3([1 size(Cmap3,1)*0.5:end],:));
    imagesc(Cmax(:,:)); colormap(Cmap3); 
    caxis([-max(Cmax_tab(~isinf(Cmax_tab)))*1.1 max(Cmax_tab(~isinf(Cmax_tab)))*1.1]);
%     FxEIT_PlotBoundary(background_mask,'w',2);
    FxEIT_PlotBoundary(ROI,'w',1);
    set(gca,'Ydir','normal'); 
    axis image off; title('C_m_a_x Image');
    subplottight(4,size(Image,3)+1,size(Image,3)+3);
    colormap(Cmap3); caxis([-max(max(Cmax(isfinite(Cmax)))) max(max(Cmax(isfinite(Cmax))))]);
    h = colorbar('Location','west'); set(h,'ylim',[0,max(max(Cmax(isfinite(Cmax))))]);
    axis off;
    
    % Display Maximum Pressure mapping
    subplottight(4,size(Image,3)+1,2*size(Image,3)+3); axis off;
    imagesc(Pstar(:,:)); colormap(Cmap3); caxis([-max(PEEP) max(PEEP)]);
%     FxEIT_PlotBoundary(background_mask,'w',2);
    FxEIT_PlotBoundary(ROI,'w',1);
    set(gca,'Ydir','normal'); 
    axis image off; title('P_m_a_x Image');
    subplottight(4,size(Image,3)+1,2*size(Image,3)+4);
    colormap(Cmap3); caxis([-max(PEEP) max(PEEP)]);
    h = colorbar('Location','west'); set(h,'ylim',[0,max(PEEP)]);
    axis off;
    
    % Display Overdist vs. Collapse plot
    subplottight(4,size(Image,3)+1,[size(Image,3)+5:2*(size(Image,3)+1)-1, 2*size(Image,3)+6:3*(size(Image,3)+1)-1]);
    plot(1:length(PEEP),CollapseSum,'-co','LineWidth',1.5); hold on; 
    plot(1:length(PEEP),OverdistSum,'-yo','LineWidth',1.5);  legend({'CollapseSum','OverdistSum'},'Location','NorthWest');
    ylim_COD = max([max(OverdistSum) max(CollapseSum)]);
    %     ylim_COD = round(ylim_COD/10)*10+10;
    ylim_COD = 100;
    ylim([0 ylim_COD]);
    set(gca,'XTick',1:length(PEEP))
    set(gca,'XTickLabel',PEEP)
    xlabel('PEEP','Color','w');
    ylabel('Percentage (%)','Color','w');
    title('< Overdistention vs. Collapse >','Color','w','Fontsize',15);
%     Image_tab = Image(:,:,1);

%% saparate display
    Cmap3 = FxImage_Cmap(3); % import double side colormap
    Cmap3 = Cmap3(end:-1:1,:); % reverse colormap
    Cmap3(1,:) = [1 1 1]; % background white
    figure; set(gcf,'units','normalized','outerposition',[0.3 0.2 0.58 0.63]);
    whitebg('k'); set(gcf,'color','k');
    for i = 1:size(Image,3)
        subplottight(3,size(Image,3)+1,i);
        imagesc(COD(:,:,i)); caxis([-1 1]); colormap(Cmap3); freezeColors;
%         set(gca,'Ydir','normal'); 
        axis image off;
%         FxEIT_PlotBoundary(background_mask,'w',2);
        FxEIT_PlotBoundary(ROI,'w',1);
        title(['PEEP : ',num2str(PEEP(i))],'Color','w');
        set(gca, 'XTickLabel', [],'XTick',[])
    end
    subplottight(3,size(Image,3)+1,size(Image,3)+1);
    caxis([-1 1]); colormap(Cmap3); colorbar('Ticks',[-1,1],'TickLabels',{'-1 (Collapse)','1 (Overdist)'},'Location','west'); 
    axis off; 
    
    % Display Maximum Compliance Image
    subplottight(3,size(Image,3)+1,size(Image,3)+2);
    % colormap(Cmap3([1 size(Cmap3,1)*0.5:end],:));
    imagesc(Cmax(:,:)); colormap(Cmap3); 
    caxis([-max(Cmax_tab(~isinf(Cmax_tab)))*1.1 max(Cmax_tab(~isinf(Cmax_tab)))*1.1]);
    FxEIT_PlotBoundary(background_mask,'w',2);
    FxEIT_PlotBoundary(ROI,'w',1);
%     set(gca,'Ydir','normal'); 
    axis image off; 
%     title('C_m_a_x Image');
    subplottight(3,size(Image,3)+1,size(Image,3)+3);
    colormap(Cmap3); caxis([-max(max(Cmax(isfinite(Cmax)))) max(max(Cmax(isfinite(Cmax))))]);
    h = colorbar('Location','west'); set(h,'ylim',[0,max(max(Cmax(isfinite(Cmax))))]);
    axis off;
    
    % Display Maximum Pressure mapping
    subplottight(3,size(Image,3)+1,2*size(Image,3)+3); axis off;
    temp = Pstar;
    temp(isinf(temp)) = -inf;
    imagesc(temp(:,:)); colormap(Cmap3); caxis([-max(PEEP) max(PEEP)]);
%     FxEIT_PlotBoundary(background_mask,'w',2);
    FxEIT_PlotBoundary(ROI,'w',1);
%     set(gca,'Ydir','normal'); 
    axis image off;
    Cmap4 = zeros(length(PEEP)*2,3);
    %     Cmap4(length(PEEP)*2:-1:length(PEEP)+1,:) = prism(length(PEEP));
    Cmap4(length(PEEP)+1:length(PEEP)*2,:) = prism(length(PEEP));
    colormap(Cmap4);
%     title('P_m_a_x Image');
    subplottight(3,size(Image,3)+1,2*size(Image,3)+4);
    colormap(Cmap4); caxis([-max(PEEP) max(PEEP)]);
    h = colorbar('Location','west'); set(h,'ylim',[0,max(PEEP)]);
    axis off;
    
    % Display Overdist vs. Collapse plot
    subplottight(3,size(Image,3)+1,[size(Image,3)+4:2*(size(Image,3)+1), 2*size(Image,3)+5:3*(size(Image,3)+1)]);
    hold on; plot(1:length(PEEP),CollapseSum,'-co','LineWidth',1.5); 
    plot(1:length(PEEP),OverdistSum,'-yo','LineWidth',1.5);  legend({'CollapseSum','OverdistSum'},'Location','NorthWest');
    ylim_COD = max([max(OverdistSum) max(CollapseSum)]);
    %     ylim_COD = round(ylim_COD/10)*10+10;
    ylim_COD = 50;
    ylim([0 ylim_COD]);
    set(gca,'XTick',1:length(PEEP))
    set(gca,'XTickLabel',PEEP)
    xlabel('PEEP','Color','w');
    ylabel('Percentage (%)','Color','w');
    title('< Overdistention vs. Collapse >','Color','w','Fontsize',15);
%     Image_tab = Image(:,:,1);
  
end

% COD
% Collapse
% Overdist
% Pstar
% Cmax
% CollapseSum
% OverdistSum