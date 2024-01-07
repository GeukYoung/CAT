function [optPEEP_COD, Result] = FxParam_optPEEP(Tidal,idx_PEEP,n_avg)
    if nargin < 3
        n_avg = 1;
    end
    PEEP_lv = round([Tidal(idx_PEEP).PEEP],1);
    nPEEP = length(PEEP_lv);
    
    if ~isfield(Tidal,'Im_Cdyn')
        if isfield(Tidal,'dP')
            for cnt = 1:length(Tidal)
                Tidal(cnt).Im_Cdyn = Tidal(cnt).Im_TV ./ Tidal(cnt).dP;
            end
        else
            error('Pressure data is empty');
        end
    end
    
    if ~isfield(Tidal,'Im_mask_lung')
        opt.th_lung = 0.25;
        for cnt = 1:length(Tidal)
            Tidal(cnt).Im_mask_lung = Tidal(cnt).Im_TV > nanmax(Tidal(cnt).Im_TV) * opt.th_lung;
        end
    end
    
    Im_COD = zeros(size(Tidal(idx_PEEP(1)).Im_Cdyn,1),nPEEP);
    Im_Pstar = zeros(size(Tidal(idx_PEEP(1)).Im_Cdyn,1),1);
    Im_Cmax = zeros(size(Tidal(idx_PEEP(1)).Im_Cdyn,1),1);
    
    
    idx_PEEPmask = [];
    for cnt = 0:n_avg-1
        idx_PEEPmask = [idx_PEEPmask idx_PEEP-cnt];
    end
    
    mask_total = max([Tidal(idx_PEEPmask).Im_mask_lung],[],2);
    idx_lung = find(mask_total);
    
    C_tab = zeros(length(Tidal(idx_PEEP(1)).Im_Cdyn),1);
    for cnt = 0:n_avg-1
        C_tab = C_tab + [Tidal(idx_PEEP-cnt).Im_Cdyn];
    end
    C_tab = C_tab/n_avg;
    C_tab = C_tab(idx_lung,:);
    
    % find Compliance_max & Pstar
    [Cmax, idx_Pstar] = max(C_tab,[],2);
    Pstar = PEEP_lv(idx_Pstar)';
    
    for cnt_PEEP = 1:nPEEP
        for cnt_pixel = 1:length(Cmax)
            if PEEP_lv(cnt_PEEP) < Pstar(cnt_pixel) % collapse
                Collapse_tab(cnt_pixel,cnt_PEEP) = ...
                    (Cmax(cnt_pixel) - C_tab(cnt_pixel,cnt_PEEP))/Cmax(cnt_pixel);
                Overdist_tab(cnt_pixel,cnt_PEEP) = 0;
            elseif PEEP_lv(cnt_PEEP) > Pstar(cnt_pixel) % overdist
                Collapse_tab(cnt_pixel,cnt_PEEP) = 0;
                Overdist_tab(cnt_pixel,cnt_PEEP) = ...
                    (Cmax(cnt_pixel) - C_tab(cnt_pixel,cnt_PEEP))/Cmax(cnt_pixel);
            else
                Collapse_tab(cnt_pixel,cnt_PEEP) = 0;
                Overdist_tab(cnt_pixel,cnt_PEEP) = 0;
            end
        end
        CollapseSum(cnt_PEEP) = round(sum(Collapse_tab(:,cnt_PEEP).*Cmax)/sum(Cmax)*100*10)/10;
        OverdistSum(cnt_PEEP) = round(sum(Overdist_tab(:,cnt_PEEP).*Cmax)/sum(Cmax)*100*10)/10;
    end
    Im_COD(idx_lung,:) = (Overdist_tab - Collapse_tab);
    Im_Pstar(idx_lung,:) = Pstar;
    Im_Cmax(idx_lung,:) = Cmax;
    
    if PEEP_lv(end) > PEEP_lv(1)
        intp_PEEP_lv = PEEP_lv(1):0.1:PEEP_lv(end);
    else
        intp_PEEP_lv = PEEP_lv(1):-0.1:PEEP_lv(end);
    end
    intp_CollapseSum = interp1(PEEP_lv,CollapseSum,intp_PEEP_lv,'linear');
    intp_OverdistSum = interp1(PEEP_lv,OverdistSum,intp_PEEP_lv,'linear');
    intp_hms = interp1(PEEP_lv,[Tidal(idx_PEEP).t_insp],intp_PEEP_lv,'linear');
    intp_Cdyn = interp1(PEEP_lv,[Tidal(idx_PEEP).Cdyn],intp_PEEP_lv,'linear');
%     figure;
%     plot(intp_PEEP_lv,intp_CollapseSum); hold on; plot(PEEP_lv,CollapseSum,'o');
%     plot(intp_PEEP_lv,intp_OverdistSum); hold on; plot(PEEP_lv,OverdistSum,'o');

    [~, tp_C] = max(intp_Cdyn);
    optPEEP_Cdyn = intp_PEEP_lv(tp_C);
    [~, tp_COD] = min(abs(intp_CollapseSum-intp_OverdistSum));
    optPEEP_COD = intp_PEEP_lv(tp_COD);
    
%% GI opt
    TV_tab = [Tidal(idx_PEEP).Im_Cdyn];
    TV_tab = TV_tab(idx_lung,:);
    for cnt_PEEP = 1:nPEEP
        Im_TV = TV_tab(:,cnt_PEEP);
        Im_TV(Im_TV<0) = 0;
        GI(cnt_PEEP) = (sum(abs(Im_TV-median(Im_TV))))/sum(Im_TV);
    end
    
%% CoV opt
    % ROI selection
    ROIx = sum(reshape(mask_total,64 ,64),2);
    idx_ROIx = [find(ROIx>0,1,'first')-1 find(ROIx>0,1,'last')+1];
    ROIy = sum(reshape(mask_total,64 ,64),1);
    idx_ROIy = [find(ROIy>0,1,'first')-1 find(ROIy>0,1,'last')+1];
    
    for cnt_PEEP = 1:nPEEP
        Im_TV = Tidal(idx_PEEP(cnt_PEEP)).Im_Cdyn;
        Im_TV(~mask_total) = 0;
        Im_TV(Im_TV<0) = 0;
        Im_TV = reshape(Im_TV,sqrt(length(Im_TV)),sqrt(length(Im_TV)))';
        Im_TV(:,[1:idx_ROIx(1) idx_ROIx(2):end]) = [];
        Im_TV([1:idx_ROIy(1) idx_ROIy(2):end],:) = [];
%         figure; imagesc(Im_TV);
        
        total_sum = sum(Im_TV(:));
        % CoVx
        temp = 0;
        for i = 1:size(Im_TV,2)
            temp = temp + nansum(Im_TV(:,i)) * (i/size(Im_TV,2));
        end
        CoVx(cnt_PEEP) = temp/total_sum;
        
        % CoVy
        temp = 0;
        for i = 1:size(Im_TV,1)
            temp = temp + nansum(Im_TV(i,:)) * (i/size(Im_TV,1));
        end
        CoVy(cnt_PEEP) = temp/total_sum;
        
    end
    

%% export
    Result.CollapseSum = CollapseSum;
    Result.OverdistSum = OverdistSum;
    Result.Cdyn = [Tidal(idx_PEEP).Cdyn];
    Result.GI = [Tidal(idx_PEEP).GI];
    Result.GI_opt = GI;
    Result.CoVx_opt = CoVx;
    Result.CoVy_opt = CoVy;
    Result.optPEEP = optPEEP_COD;
    Result.optPEEP_hms = intp_hms(tp_COD);
    Result.optPEEP_pt = 0.5*(intp_CollapseSum(tp_COD)+intp_OverdistSum(tp_COD));
    Result.optPEEP_Cdyn = optPEEP_Cdyn;
    Result.Im_COD = Im_COD;
    Result.Im_Pstar = Im_Pstar;
    Result.Im_Cmax = Im_Cmax;
    Result.Im_Cdyn = [Tidal(idx_PEEP).Im_Cdyn];
end
