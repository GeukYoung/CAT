function [Tidal,DataSet] = FxParam_AirTom_v1_peakonly(DataSet)
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
        
        Tidal(cnt).TVi = sum(abs(flow_temp(flow_temp>0)))/60/100;
        Tidal(cnt).TVe = sum(abs(flow_temp(flow_temp<0)))/60/100;
        Tidal(cnt).MV = Tidal(cnt).TVi * Tidal(cnt).RR;
        Tidal(cnt).dFRC = Tidal(cnt).TVi - Tidal(cnt).TVe;
        Tidal(cnt).Cdyn = Tidal(cnt).TVi/Tidal(cnt).dP;
        Tidal(cnt).AU2mL = Tidal(cnt).TVi/Tidal(cnt).TVi_au;        
        Tidal(cnt).unit_V = 'mL';
        clear flow_temp;
    end
end

% for cnt = 1:length(peaks.insp)-1
%     %% TV
%     if 1
%         Tidal(cnt).Im_TV = DataSet.EIT.Image.RM*-(DataSet.EIT.ohm_HPF(:,peaks.insp(cnt))-DataSet.EIT.ohm_HPF(:,peaks.exp(cnt)))*median([Tidal.AU2mL]);
%         Tidal(cnt).TVi_au2ml_same = Tidal(cnt).TVi_au*median([Tidal.AU2mL]);
%     end
%     if 0
%         Tidal(cnt).Im_TV = DataSet.EIT.Image.RM*-(DataSet.EIT.ohm_HPF(:,peaks.insp(cnt))-DataSet.EIT.ohm_HPF(:,peaks.exp(cnt)))*Tidal(cnt).AU2mL;
%     end
%     Tidal(cnt).Im_mask_lung = Tidal(cnt).Im_TV > nanmax(Tidal(cnt).Im_TV) * opt.th_lung;
%     
%     %% Cdyn
%     if isfield(DataSet,'Pressure')
%         Tidal(cnt).Im_Cdyn = Tidal(cnt).Im_TV ./ Tidal(cnt).dP;
%     end
%     
%     %% RVD
%     tp = [Tidal(cnt).idx_exp Tidal(cnt).idx_insp];
%     [Tidal(cnt).Im_RVD, Tidal(cnt).sdRVD] = sFxEIT_RVD(DataSet,tp);
%     
%     %% Popen
%     if isfield(DataSet,'Pressure')
%         tp = [Tidal(cnt).idx_exp Tidal(cnt).idx_insp];
%         [Tidal(cnt).Im_Popen, Tidal(cnt).sdPopen] = sFxEIT_Popen(DataSet,tp);
%     end
%     
%     %% GI
%     Tidal(cnt).GI =  sFxEIT_GI(Tidal(cnt).Im_TV, DataSet.EIT.Image.isbnd);
%     
%     %% CoVx,y
%     [Tidal(cnt).CoVx, Tidal(cnt).CoVy] = sFxEIT_CoV(Tidal(cnt).Im_TV);
%     
%     %% avg filter
% %     if cnt > opt.n_avg-1
% %         Tidal(cnt).sdRVD = nanmean([Tidal(cnt-opt.n_avg+1:cnt).sdRVD]);
% %         Tidal(cnt).sdPopen = nanmean([Tidal(cnt-opt.n_avg+1:cnt).sdPopen]);
% %         Tidal(cnt).GI = nanmean([Tidal(cnt-opt.n_avg+1:cnt).GI]);
% %         Tidal(cnt).CoVx = nanmean([Tidal(cnt-opt.n_avg+1:cnt).CoVx]);
% %         Tidal(cnt).CoVy = nanmean([Tidal(cnt-opt.n_avg+1:cnt).CoVy]);
% %         Tidal(cnt).MV = nanmean([Tidal(cnt-opt.n_avg+1:cnt).MV]);
% %     end
% end

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

    function [Im_RVD, sdRVD] = sFxEIT_RVD(DataSet,tp)
        opt.th_lung = 0.25;
        opt.th_rvd = 0.4;
        
        Im_RVD = zeros(size(DataSet.EIT.Image.RM,1),1);
        temp_TV = DataSet.EIT.Image.RM*-(DataSet.EIT.ohm_HPF(:,tp(2))-DataSet.EIT.ohm_HPF(:,tp(1)));
        Im_mask_lung = temp_TV > (nanmax(temp_TV)*opt.th_lung);
        idx_lung = find(Im_mask_lung);
        wave_sigma = DataSet.EIT.Image.RM(idx_lung,:)*-(DataSet.EIT.ohm_HPF(:,tp(1):tp(2))-DataSet.EIT.ohm_HPF(:,tp(1)));

        wave_sigma_sum = sum(wave_sigma);
        t_global = find(wave_sigma_sum>(max(wave_sigma_sum)*opt.th_rvd),1,'first')/DataSet.EIT.fs*1000; % s -> ms;
        for cnt_RVD = 1:length(idx_lung)
            t_open(cnt_RVD) = find(wave_sigma(cnt_RVD,:)>(max(wave_sigma(cnt_RVD,:))*opt.th_rvd),1,'first')/DataSet.EIT.fs*1000; % s -> ms
            Im_RVD(idx_lung(cnt_RVD)) = t_open(cnt_RVD) - t_global;
        end
        sdRVD = std(t_open);
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