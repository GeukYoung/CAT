function [Tidal,DataSet] = FxParam_AirTom_v1(DataSet)
% input
%   data : EIT data 
% output
%   tidal : tidal volume & time info

DataSet.EIT.fs = round(1/seconds(median(diff(DataSet.EIT.t_hms))));
opt.th_lung = 0.2;
opt.n_avg = 3;
%% HPF
if ~isfield(DataSet.EIT,'ohm_HPF')
    Wn = 1; N = 2;
    [b,a] = butter(N, Wn*2/DataSet.EIT.fs,'low');
    
    idx_notfin = ~isfinite(DataSet.EIT.ohm);
    DataSet.EIT.ohm_HPF = fillmissing(DataSet.EIT.ohm','nearest','EndValues','nearest')';
    for cnt = 1:size(DataSet.EIT.ohm)
        DataSet.EIT.ohm_HPF(cnt,:) = filtfilt(b,a,DataSet.EIT.ohm_HPF(cnt,:));
    end
    DataSet.EIT.ohm_HPF(idx_notfin) = NaN;
    clear a b Wn N idx_notfin;
end

%% RVS
if ~isfield(DataSet.EIT,'RVS')
    DataSet.EIT.RVS = -sum(DataSet.EIT.Image.RM) * DataSet.EIT.ohm_HPF;
end

%% Check Saturation
opt.th_sat = 15;
opt.w_sat = 2*100;
DataSet.EIT.Sat = (mean(DataSet.EIT.ohm_HPF) > opt.th_sat) | isnan(mean(DataSet.EIT.ohm_HPF));
DataSet.EIT.Sat = movmax(DataSet.EIT.Sat,[opt.w_sat opt.w_sat]);
DataSet.EIT.RVS(DataSet.EIT.Sat) = nan;

%% Find peaks
peaks = sFxEIT_TVpeak(DataSet.EIT.RVS,DataSet.EIT.fs);
DataSet.EIT.peaks = peaks;

for cnt = 1:length(peaks.insp)-1
    Tidal(cnt).idx_exp = peaks.exp(cnt);
    Tidal(cnt).idx_insp = peaks.insp(cnt);
%     Tidal(cnt).idx_exp_v = peaks.exp(cnt):(peaks.insp(cnt+1)-1);
%     Tidal(cnt).idx_insp_v = peaks.insp(cnt):(peaks.exp(cnt)-1);
    Tidal(cnt).t_exp = DataSet.EIT.t_hms(peaks.exp(cnt));
    Tidal(cnt).t_insp = DataSet.EIT.t_hms(peaks.insp(cnt));
    Tidal(cnt).td_exp = DataSet.EIT.t_hms(peaks.exp(cnt+1)) - DataSet.EIT.t_hms(peaks.insp(cnt));
    Tidal(cnt).td_insp = DataSet.EIT.t_hms(peaks.insp(cnt)) - DataSet.EIT.t_hms(peaks.exp(cnt));
    Tidal(cnt).td_tidal = Tidal(cnt).td_insp + Tidal(cnt).td_exp;
    Tidal(cnt).RR = 60/seconds(Tidal(cnt).td_tidal);
    Tidal(cnt).IE = seconds(Tidal(cnt).td_insp)/seconds(Tidal(cnt).td_tidal);
    Tidal(cnt).TVi_au = DataSet.EIT.RVS(peaks.insp(cnt)) - DataSet.EIT.RVS(peaks.exp(cnt));
    Tidal(cnt).TVe_au = DataSet.EIT.RVS(peaks.insp(cnt)) - DataSet.EIT.RVS(peaks.exp(cnt+1));
    Tidal(cnt).MV = Tidal(cnt).TVi_au * Tidal(cnt).RR;
    Tidal(cnt).dFRC = Tidal(cnt).TVi_au - Tidal(cnt).TVe_au;
    Tidal(cnt).AU2mL = 1;
    Tidal(cnt).unit_V = 'AU';
    
    if isfield(DataSet,'Pressure')
        Tidal(cnt).PIP = max(DataSet.Pressure.Paw(peaks.exp(cnt):peaks.exp(cnt+1)));
        Tidal(cnt).PEEP = min(DataSet.Pressure.Paw(peaks.exp(cnt):peaks.exp(cnt+1)));
        Tidal(cnt).dP = Tidal(cnt).PIP - Tidal(cnt).PEEP;
        
        flow_temp = DataSet.Pressure.Faw(peaks.exp(cnt):peaks.exp(cnt+1));
        Tidal(cnt).flow_offset = mean(flow_temp([1:5 end-4:end]));
%         flow_temp = flow_temp - Tidal(cnt).flow_offset;
        
        Tidal(cnt).TVi = sum(abs(flow_temp(flow_temp>0)))/60/100*1000;
        Tidal(cnt).TVe = sum(abs(flow_temp(flow_temp<0)))/60/100*1000;
        Tidal(cnt).MV = Tidal(cnt).TVi * Tidal(cnt).RR;
        Tidal(cnt).dFRC = Tidal(cnt).TVi - Tidal(cnt).TVe;
        Tidal(cnt).Cdyn = Tidal(cnt).TVi/Tidal(cnt).dP;
        Tidal(cnt).AU2mL = Tidal(cnt).TVi/Tidal(cnt).TVi_au;        
        Tidal(cnt).unit_V = 'mL';
        clear flow_temp;
    end
end

for cnt = 1:length(peaks.insp)-1
    %% TV
    if 1
        Tidal(cnt).Im_TV = DataSet.EIT.Image.RM*-(DataSet.EIT.ohm_HPF(:,peaks.insp(cnt))-DataSet.EIT.ohm_HPF(:,peaks.exp(cnt)))*median([Tidal.AU2mL]);
        Tidal(cnt).TVi_au2ml_same = Tidal(cnt).TVi_au*median([Tidal.AU2mL]);
    end
    if 0
        Tidal(cnt).Im_TV = DataSet.EIT.Image.RM*-(DataSet.EIT.ohm_HPF(:,peaks.insp(cnt))-DataSet.EIT.ohm_HPF(:,peaks.exp(cnt)))*Tidal(cnt).AU2mL;
    end
    Tidal(cnt).Im_mask_lung = Tidal(cnt).Im_TV > nanmax(Tidal(cnt).Im_TV) * opt.th_lung;
   
    %% RVD
    tp = [peaks.exp(cnt) peaks.insp(cnt)];
    [Tidal(cnt).Im_RVD, Tidal(cnt).sdRVD, Tidal(cnt).mRVD] = sFxEIT_RVD(DataSet,tp);
    
    %% Pendelluft volume
    tp = [peaks.exp(cnt) peaks.insp(cnt) peaks.exp(cnt+1)];
    [Tidal(cnt).Im_Pendelluft, Tidal(cnt).Vpendelluft, Tidal(cnt).Im_tPendelluft] = sFxEIT_Pendelluft(DataSet,tp);
    if 1
        Tidal(cnt).Im_Pendelluft = Tidal(cnt).Im_Pendelluft.*median([Tidal.AU2mL]);
        Tidal(cnt).Vpendelluft = Tidal(cnt).Vpendelluft.*Tidal(cnt).AU2mL;
    end
    
    %% Cdyn
    if isfield(DataSet,'Pressure')
        Tidal(cnt).Im_Cdyn = Tidal(cnt).Im_TV ./ Tidal(cnt).dP;
    end
    
    %% Popen
    if isfield(DataSet,'Pressure')
        tp = [peaks.exp(cnt) peaks.insp(cnt)];
        [Tidal(cnt).Im_Popen, Tidal(cnt).sdPopen] = sFxEIT_Popen(DataSet,tp);
    end
    
    %% GI
    Tidal(cnt).GI =  sFxEIT_GI(Tidal(cnt).Im_TV, DataSet.EIT.Image.isbnd);
    
    %% CoVx,y
    [Tidal(cnt).CoVx, Tidal(cnt).CoVy] = sFxEIT_CoV(Tidal(cnt).Im_TV);
    
    %% avg filter
    if cnt > opt.n_avg-1
        Tidal(cnt).sdRVD = nanmean([Tidal(cnt-opt.n_avg+1:cnt).sdRVD]);
        Tidal(cnt).sdPopen = nanmean([Tidal(cnt-opt.n_avg+1:cnt).sdPopen]);
        Tidal(cnt).GI = nanmean([Tidal(cnt-opt.n_avg+1:cnt).GI]);
        Tidal(cnt).CoVx = nanmean([Tidal(cnt-opt.n_avg+1:cnt).CoVx]);
        Tidal(cnt).CoVy = nanmean([Tidal(cnt-opt.n_avg+1:cnt).CoVy]);
        Tidal(cnt).MV = nanmean([Tidal(cnt-opt.n_avg+1:cnt).MV]);
    end
end

%% Sub functions
    function [CoVx, CoVy] = sFxEIT_CoV(Im_TV)
        Im_TV(Im_TV<0) = 0;
        total_sum = sum(Im_TV);
        Im_TV = reshape(Im_TV,sqrt(length(Im_TV)),sqrt(length(Im_TV)))';
        
        % CoVx
        temp = 0;
        for i = 1:size(Im_TV,2)
            temp = temp + nansum(Im_TV(:,i)) * (i/size(Im_TV,2));
        end
        CoVx = temp/total_sum;
        
        % CoVy
        temp = 0;
        for i = 1:size(Im_TV,1)
            temp = temp + nansum(Im_TV(i,:)) * (i/size(Im_TV,1));
        end
        CoVy = temp/total_sum;
    end

    function [GI] = sFxEIT_GI(Im_TV, isbnd)
        Im_TV(Im_TV<0) = 0;
        Im_TV(isbnd) = [];
        GI = (sum(abs(Im_TV-median(Im_TV))))/sum(Im_TV); 
    end

    function [Im_RVD, RVDI, mRVD] = sFxEIT_RVD(DataSet,tp)
        opt.th_lung = 0.25;
        opt.th_rvd = 0.4;
        
        Im_RVD = zeros(size(DataSet.EIT.Image.RM,1),1);
        temp_TV = DataSet.EIT.Image.RM*-(DataSet.EIT.ohm_HPF(:,tp(2))-DataSet.EIT.ohm_HPF(:,tp(1)));
        Im_mask_lung = temp_TV > (nanmax(temp_TV)*opt.th_lung);
        idx_lung = find(Im_mask_lung);
        wave_sigma = DataSet.EIT.Image.RM(idx_lung,:)*-(DataSet.EIT.ohm_HPF(:,tp(1):tp(2))-DataSet.EIT.ohm_HPF(:,tp(1)));

%         wave_sigma_sum = sum(wave_sigma);
%         RVD_global = find(wave_sigma_sum>(max(wave_sigma_sum)*opt.th_rvd),1,'first')/length(wave_sigma_sum); % percentage
        for cnt_RVD = 1:length(idx_lung)
            Im_RVD(idx_lung(cnt_RVD)) = find(wave_sigma(cnt_RVD,:)>(max(wave_sigma(cnt_RVD,:))*opt.th_rvd),1,'first')/length(wave_sigma(cnt_RVD,:)); % percentage
        end
        RVDI = std(Im_RVD(idx_lung));
        mRVD = mean(Im_RVD(idx_lung));
%         figure;
%         imagesc(reshape(Im_RVD,64,64)'); axis image off;
    end

    function [Im_Pendelluft, Vpendelluft, Im_tPendelluft] = sFxEIT_Pendelluft(DataSet,tp)
        opt.th_lung = 0.25;
        
        Im_tPendelluft = zeros(size(DataSet.EIT.Image.RM,1),1);
        Im_Pendelluft = zeros(size(DataSet.EIT.Image.RM,1),1);
        temp_TV = DataSet.EIT.Image.RM*-(DataSet.EIT.ohm_HPF(:,tp(2))-DataSet.EIT.ohm_HPF(:,tp(1)));
        Im_mask_lung = temp_TV > (nanmax(temp_TV)*opt.th_lung);
        idx_lung = find(Im_mask_lung); % lung ROI
        wave_sigma = DataSet.EIT.Image.RM(idx_lung,:)*-(DataSet.EIT.ohm_HPF(:,tp(1):tp(3))-DataSet.EIT.ohm_HPF(:,tp(1)));
        [~, idx_sExpG] = max(sum(wave_sigma)); % global peak time
        for cnt_Pen = 1:length(idx_lung)
            [~, idx_sExpR] = max(wave_sigma(cnt_Pen,:));
            Im_Pendelluft(idx_lung(cnt_Pen)) = wave_sigma(cnt_Pen,idx_sExpR) - wave_sigma(cnt_Pen,idx_sExpG);
            Im_tPendelluft(idx_lung(cnt_Pen)) = (idx_sExpR - idx_sExpG)/DataSet.EIT.fs;
            if Im_tPendelluft(idx_lung(cnt_Pen)) > 0 % Gas Lost (after global peak)
                Im_Pendelluft(idx_lung(cnt_Pen)) =  -Im_Pendelluft(idx_lung(cnt_Pen)); % sign code
            end
        end
        Vpendelluft = sum(abs(Im_Pendelluft(idx_lung)));
        % positive: Gained, negative: Lost
    end

    function [Im_Popen, sdPopen] = sFxEIT_Popen(DataSet,tp)
        opt.th_lung = 0.25;
        opt.th_popen = 0.1;
        
        Im_Popen = zeros(size(DataSet.EIT.Image.RM,1),1);
        temp_TV = DataSet.EIT.Image.RM*-(DataSet.EIT.ohm_HPF(:,tp(2))-DataSet.EIT.ohm_HPF(:,tp(1)));
        Im_mask_lung = temp_TV > (nanmax(temp_TV)*opt.th_lung);
        idx_lung = find(Im_mask_lung);
        wave_sigma = DataSet.EIT.Image.RM(idx_lung,:)*-(DataSet.EIT.ohm_HPF(:,tp(1):tp(2))-DataSet.EIT.ohm_HPF(:,tp(1)));
        for cnt_RVD = 1:length(idx_lung)
            idx_open(cnt_RVD) = find(wave_sigma(cnt_RVD,:)>(max(wave_sigma(cnt_RVD,:))*opt.th_popen),1,'first')-1;
            Im_Popen(idx_lung(cnt_RVD)) = DataSet.Pressure.Paw(tp(1) + idx_open(cnt_RVD));
        end
        sdPopen = std(Im_Popen(idx_lung));
        
%         figure;
%         imagesc(reshape(Im_Popen,64,64)'); axis image off;
    end
            
    function peaks = sFxEIT_TVpeak(RVS,fs)
        ws = round(1*fs);
        MAC = zeros(length(RVS),1);
        MAC(1:ws) = mean(RVS(1:ws));
        if MAC(ws) < RVS(ws)
            state_itc = 1;
        else
            state_itc = -1;
        end
        margin = 0.1;
        cnt_itc = ws;
        cnt_peak = 1;
        cnt_v = 1;
        i_itc = 0;
        roc = 0;
        
        for i = (ws+1):length(RVS)
            if i < 5*ws + 1
                MAC(i) = mean(RVS((i-ws):i)) - state_itc*std(RVS((i-ws):i))*margin;
            else
                MAC(i) = mean(RVS((i-ws):i)) - state_itc*std(RVS((i-5*ws):i))*margin;
            end
            
            if state_itc == 1
                if MAC(i) >= RVS(i)
                    [~,temp_roc] = max(RVS((i-cnt_itc):i));
                    roc(cnt_peak) = i - cnt_itc + temp_roc - 1;
                    i_itc(cnt_peak) = i;
                    cnt_peak = cnt_peak + 1;
                    cnt_itc = 0;
                    state_itc = -1;
                end
            elseif state_itc == -1
                if MAC(i) <= RVS(i)
                    [~,temp_roc] = min(RVS((i-cnt_itc):i));
                    roc(cnt_peak) = i - cnt_itc + temp_roc - 1;
                    i_itc(cnt_peak) = i;
                    cnt_peak = cnt_peak + 1;
                    cnt_itc = 0;
                    state_itc = 1;
                end
            end
            cnt_itc = cnt_itc + 1;
        end
        
        locs_Peak2 = roc;
        if RVS(locs_Peak2(1)) > RVS(locs_Peak2(2))
            locs_Peak2(1) = [];
        end
        if mod(length(locs_Peak2),2) == 1
            locs_Peak2(end) = [];
        end
        
        tidal_th = 0.3;
        idx_rm = [];
        for i = 2:(length(locs_Peak2)/2)
            pre_tidal = abs(RVS(locs_Peak2(2*(i-1)-1)) - RVS(locs_Peak2(2*(i-1))));
            cur_tidal = abs(RVS(locs_Peak2(2*i-1)) - RVS(locs_Peak2(2*i)));
            if (tidal_th*pre_tidal) > cur_tidal
                idx_rm = [idx_rm i];
            end
        end
        
        if length(idx_rm) > 1
            locs_Peak2([(2*idx_rm-1) (2*idx_rm)]) = [];
        end
        peaks.exp = locs_Peak2(1:2:end)';
        peaks.insp = locs_Peak2(2:2:end)';
    end
end
%% test image
% figure; plot(DataSet.EIT.RVS)
% tp = [119432 119515];
% 
% temp = -(DataSet.EIT.ohm_HPF(:,tp(2)) - DataSet.EIT.ohm_HPF(:,tp(1)));
% temp2 = DataSet.EIT.Image.RM * temp;
% temp2 = reshape(temp2,64,64)';
% 
% F = griddedInterpolant({1:64 1:64},temp2);
% 
% imagesc(F({linspace(1,64,2048) linspace(1,64,2048)})); colormap(Image_AirTom.Cmap_oneside);
% hold on; imagesc(Image_AirTom.mask_cdata,'AlphaData',Image_AirTom.mask_alpha);
% axis image off;
% 
% Image_AirTom.Cmap_oneside = Cmap_oneside;
% Image_AirTom.Cmap_doubleside = Cmap_doubleside;
% Image_AirTom.Im_size = 2048;


%%

% function [tidal] = FxEIT_findRespPeak(data,time)
% % input
% %   data : EIT data 
% %   time : EIT time vector
% %   PSG : PSG data
% % output
% %   tidal : tidal volume & time info
% 
% if nargin < 2
%     time = 1:length(data);
% end
% 
% fs = 1/mean(diff(time)); % EIT framerate
% if size(data,1) == 1
%     data = data';
% end
% 
% %% detrend data
% if nargin < 2
%     Detrend_data = data; % filtering
% else
%     N = 1; Wn = 0.8;
%     Detrend_data = FxEIT_Filter(data',fs,N,Wn,1,1)'; % filtering
% end
% dd_data = diff(diff(Detrend_data)); % to find curve shape
% 
% %% find insp/exp peak
% [~,locs_temp1] = findpeaks(Detrend_data,...
%     'MinPeakDistance',round(fs*1)); % find vally
% [~,locs_temp2] = findpeaks(-Detrend_data,...
%     'MinPeakDistance',round(fs*1)); % find peak
% locs_Peak = sort([locs_temp1' locs_temp2']); % combine both peak data
% clear locs_temp1 locs_temp2;
% 
% amp_Peak = data(locs_Peak); % amplitude at peak point
% sign_dd_data = sign(dd_data(locs_Peak)); % to find curve shape
% [a, ~] = find(sign_dd_data==1); % find first exp point
% 
% if a(1) > 1 % exp is start
%     locs_Peak(1:a-1) = [];
%     amp_Peak(1:a-1) = [];
%     sign_dd_data(1:a-1) = [];
% end
% 
% % figure; plot(Detrend_data);
% % hold on; plot(locs_Peak(1:2:end),Detrend_data(locs_Peak(1:2:end)),'^b', ...
% %      'MarkerFaceColor',[0 0 1],'MarkerSize',5);
% % hold on; plot(locs_Peak(2:2:end),Detrend_data(locs_Peak(2:2:end)),'vr', ...
% %      'MarkerFaceColor',[1 0 0],'MarkerSize',5);
% % xlim([5500 7500]); ylim([9200000 10300000]);
% % 
% % figure; plot(Detrend_data);
% % hold on; plot(locs_Peak(1:2:end),Detrend_data(locs_Peak(1:2:end)),'vr', ...
% %      'MarkerFaceColor',[1 0 0],'MarkerSize',5);
% % hold on; plot(locs_Peak(2:2:end),Detrend_data(locs_Peak(2:2:end)),'^b', ...
% %      'MarkerFaceColor',[0 0 1],'MarkerSize',5);
% % xlim([0 1000500]); ylim([9200000 10.300000]);
% 
% %% insp/exp extraction 
% cnt_idx = 1;
% sign_curve = 1;
% for i = 1:length(locs_Peak)-1
%     if sign_dd_data(i) ~= sign_curve
%         if sign_dd_data(i) > 0
%             if amp_Peak(i) < amp_Peak(idx_locs(cnt_idx-1))
%                 idx_locs(cnt_idx-1) = i;
%             end
%         else sign_dd_data(i) < 0;
%             if amp_Peak(i) > amp_Peak(idx_locs(cnt_idx-1))
%                 idx_locs(cnt_idx-1) = i;
%             end
%         end
%     else
%         idx_locs(cnt_idx) = i;
%         cnt_idx = cnt_idx + 1;
%         sign_curve = -sign_curve;
%     end
% end
% locs_Peak = locs_Peak(idx_locs);
% amp_Peak = amp_Peak(idx_locs);
% sign_dd_data = sign_dd_data(idx_locs);
% clear cnt_idx sign_curve idx_locs;
% 
% if mod(length(locs_Peak),2) ~= 0 % remove last exp point
%     locs_Peak(end) = [];
%     amp_Peak(end) = [];
%     sign_dd_data(end) = [];
% end
% 
% % fig = figure; set(fig,'units','normalized' ,'outerposition' ,[0.15 0.20 0.36 0.31]); 
% % for i = 1:(length(locs_Peak)/2)-1
% %     hold on; plot(([locs_Peak(2*i-1):locs_Peak(2*i)])/fs ...
% %         ,data(locs_Peak(2*i-1):locs_Peak(2*i)),'k', 'LineWidth',2);
% %    hold on; plot(([locs_Peak(2*i):locs_Peak(2*i+1)])/fs ...
% %         ,data(locs_Peak(2*i):locs_Peak(2*i+1)),'Color',[0.2 0.9 0.5], 'LineWidth',2);
% % end
% % hold on; plot((locs_Peak(1:2:end))/fs,data(locs_Peak(1:2:end)),'b^', ...
% %      'MarkerFaceColor',[0 0 1],'MarkerSize',5);
% % hold on; plot((locs_Peak(2:2:end))/fs,data(locs_Peak(2:2:end)),'rv', ...
% %      'MarkerFaceColor',[1 0 0],'MarkerSize',5);
% % set(gca,'ytick',[]); xlabel('Time(s)','fontsize',8); ylabel('Volume (AU)','fontsize',10);
% % xlim([7501 9000]);
% 
% tidal.idx_exp = locs_Peak(1:2:end)';
% tidal.idx_insp = locs_Peak(2:2:end)';
% tidal.t_exp = time(tidal.idx_exp);
% tidal.t_insp = time(tidal.idx_insp);
% tidal.volume_insp = abs(data(tidal.idx_insp) - data(tidal.idx_exp));
% tidal.volume_exp = abs(data(tidal.idx_insp(1:end-1)) - data(tidal.idx_exp(2:end)));
% tidal.volume_residual = tidal.volume_insp(1:end-1) - tidal.volume_exp;
% tidal.interval = [0; (time(tidal.idx_insp(2:end))-time(tidal.idx_insp(1:end-1)))];
% tidal.td_insp = time(tidal.idx_insp) - time(tidal.idx_exp);
% tidal.td_exp = time(tidal.idx_exp(2:end)) - time(tidal.idx_insp(1:end-1));
% tidal.td_tidal = tidal.td_insp(1:end-1) + tidal.td_exp;
% tidal.IE = tidal.td_insp(1:end-1)./tidal.td_tidal;
% tidal.flow_insp = tidal.volume_insp./tidal.td_insp;
% tidal.flow_exp = tidal.volume_exp./tidal.td_exp;
% tidal.RR = 60./(time(tidal.idx_insp(2:end))-time(tidal.idx_insp(1:end-1)));
% 
% fig = figure; set(fig,'units','normalized' ,'outerposition' ,[0.15 0.20 0.4 0.36]);
% subplot(121); histogram(tidal.volume_insp,20,'BinWidth',std(tidal.volume_insp)/5.5,'facecolor',[0.85 0.33 0.1]);
% title('Tidal Volume'); xlabel('Volume(AU)'); ylabel('Number of tidal');
% % hold on; stem(Tidalvolume_m,20,'r');
% subplot(122); histogram(tidal.RR,20,'BinWidth',std(tidal.RR)/3.5,'facecolor',[0.85 0.33 0.1]); 
% title('Respiration Rate'); xlabel('Rate(/min)'); ylabel('Number of tidal');
% % hold on; stem(Tidalinterval_m,20,'r');
% 
% figure; set(gcf,'units','normalized' ,'outerposition' ,[0 0.20 1 0.6]);
% fig2{1} = subplottight(2,1,1,2);
% plot(time,Detrend_data);
% hold on; plot(time(locs_Peak(1:2:end)),Detrend_data(locs_Peak(1:2:end)),'^b', ...
%      'MarkerFaceColor',[0 0 1],'MarkerSize',5);
% hold on; plot(time(locs_Peak(2:2:end)),Detrend_data(locs_Peak(2:2:end)),'vr', ...
%      'MarkerFaceColor',[1 0 0],'MarkerSize',5);
% xlim([-inf inf]); set(gca,'xtick',[]); set(gca,'ytick',[]);
% ylabel('EIT Pixelsum','fontsize',15,'Rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','right');
% 
% fig2{2} = subplottight(2,1,2,2);
% plot(tidal.t_insp,tidal.volume_insp,'-o','linewidth',1.5,'markersize',2); 
% xlim([-inf inf]); set(gca,'xtick',[]); set(gca,'ytick',[]);
% ylabel('Tidal Volume','fontsize',15,'Rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','right');
% linkaxes([fig2{:}],'x');