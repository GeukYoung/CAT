% axisAll = gca;
colorMat =  [ 1 0 0];
colorMat2 = [ 0 0 1];
colorMat3 = [ 0 1 0];

%
AlphaTrans = 0.2;

textRatio = 0.1;


% for axTmp = [ ax1 ax2 ax3 ax4 ax5]
cnt = 1;
for axTmp = axisAll
    
%     if isempty(axTmp.YLabel.String)
%         continue;
%     end
    
    if isempty(axTmp)
        continue;
    end
    
    
    ylimTmp = axTmp.YLim;
    
    y1 = ylimTmp(1);
    y2 = ylimTmp(2);
    
    
    for k = 1:size(intervention,1)
        x1 = intervention{k,1}(1);
        x2 = intervention{k,1}(2);
        
        str = intervention{k,2};
        
        hfill=fill(axTmp,[x1 x2 x2 x1],[y1 y1 y2 y2 ],intervention{k,3}, 'FaceAlpha',intervention{k,4},'EdgeColor','none');
        if isShowLabel && cnt==1
            text(axTmp, x1,y1*0.2+y2*0.8,str,'fontsize',15)
        end
        uistack(hfill,'bottom')
    end
    
    cnt=cnt+1;
end
