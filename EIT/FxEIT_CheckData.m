function [temp] = FxEIT_CheckData(Data,imdl)
if nargin < 2
    imdl = FxEIT_FER2;
end
op.fs = 100;
op.ch = 16;
op.CI_max = 10000;
op.CI_ch=[2:17:256 241];
op.S_idx = 1;

disp.CI = 1;

%% 1) calculate CI
if any(Data.CI_flag)
    CI_flag = Data.CI_flag;
    [~, CI_idx]=find(CI_flag==1);

    % unstable data detection (gain change)
    %     op.S_idx=find(Data.EIT1.AS_flag==max(Data.EIT1.AS_flag),1,'first');
    op.S_idx = 1;
    if op.S_idx<(max(CI_idx)+1)
        op.S_idx=max(CI_idx);
    end

    CI=Data.EIT_v(op.CI_ch,CI_idx)./Data.Amp(CI_idx); % contact impedance
    CI2 = CI;
    CI2(repmat(sum(CI>op.CI_max)>1,16,1)) = NaN;
    CI2(:,logical(sum(isnan(CI2)))) = [];

    result.CI_raw = CI2;
    if size(result.CI_raw,2) > 1
        result.CI_min = min(result.CI_raw');
        result.CI_max = max(result.CI_raw');
        result.CI_avg = mean(result.CI_raw');
        result.CI_std = std(result.CI_raw');
    else
        result.CI_min = result.CI_raw;
        result.CI_max = result.CI_raw;
        result.CI_avg = result.CI_raw;
        result.CI_std = 0;
    end
else
    disp.CI = 0;
end
clear CI CI2;

%% 2) calculate RE
RE=FxEIT_RE(Data.EIT_v,FxEIT_mask(16));
result.RE_raw = RE.RE_raw;
result.RE_min = nanmean(RE.RE_min);
result.RE_max = nanmean(RE.RE_max);
result.RE_avg = nanmean(RE.RE_avg);
op.idx_full = 1:op.ch^2;
op.idx_full(FxEIT_mask(16)) = [];
temp_full(op.idx_full) = nanmean(RE.RE_raw');
temp_full(FxEIT_mask(op.ch)) = NaN;
temp_full = reshape(temp_full,op.ch,op.ch);
result.RE_16ch = nanmean(temp_full);
clear RE temp_full;

%% 3) post result
Data.EIT_v(:,1:op.S_idx+500) = [];
Data.Amp(:,1:op.S_idx+500) = [];
Data.Gain(:,1:op.S_idx+500) = [];

% ) ohm scale & pixel sum
Data.EIT_v_ohm = Data.EIT_v./repmat(Data.Amp,256,1);
Data.EIT_v_ohm(FxEIT_mask(16),:) = [];
Data.EIT_v_ohm_mean = mean(Data.EIT_v_ohm);
Data.pixel_sum = -ones(1,size(imdl.solve_use_matrix.RM,1)) * imdl.solve_use_matrix.RM * imdl.Proj_Mat * Data.EIT_v_ohm;
% Data.pixel_sum = -ones(1,size(imdl.solve_use_matrix.RM,1)) * imdl.solve_use_matrix.RM * Data.EIT_v_ohm;

N = 8; fc = 1; butter = 1; bpf = 1;% guk fc=8 % fc=0.5 previously
Data.pixel_sum_fc1hz = FxEIT_Filter(Data.pixel_sum,op.fs,N,fc,butter,bpf);
clear N fc butter bpf;
% figure; plot(Data.EIT1.pixel_sum); hold on; plot(Data.EIT1.pixel_sum_fc1hz,'r');

% ) find stable section
TV1 = FxEIT_findRespPeak3(Data.pixel_sum_fc1hz,(1:length(Data.pixel_sum_fc1hz))/op.fs);
TV1.RR_full = nan(length(Data.pixel_sum_fc1hz),1);
TV1.RR_full(TV1.idx_insp(2:end)) = TV1.RR;
TV1.RR_full = fillmissing(TV1.RR_full,'next');
TV1.TV_full = nan(length(Data.pixel_sum_fc1hz),1);
TV1.TV_full(TV1.idx_insp) = TV1.volume_insp;
TV1.TV_full = fillmissing(TV1.TV_full,'next');

op.find_ws = 6000;
cnt = 1;
clear TV1.std_TV TV1.std_RR;
for i = 1:op.fs:length(Data.pixel_sum_fc1hz)-op.find_ws*10
    if any(isnan(TV1.TV_full(i:(i+op.find_ws-1))))
        TV1.std_TV(cnt) = NaN;
    else
        TV1.std_TV(cnt) = std(TV1.TV_full(i:(i+op.find_ws-1)));
    end
    
    if any(isnan(TV1.RR_full(i:(i+op.find_ws-1))))
        TV1.std_RR(cnt) = NaN;
    else
        if any(TV1.RR_full(i:(i+op.find_ws*1.1)) < 8)
            TV1.std_RR(cnt) = NaN;
        else
            TV1.std_RR(cnt) = nanstd(TV1.RR_full(i:(i+op.find_ws-1)));
        end
    end
    cnt = cnt + 1;
end

% ) find Zinsp Zexp
result.Z_exp = max(mean(Data.EIT_v_ohm(:,TV1.idx_exp)'));
result.Z_insp = max(mean(Data.EIT_v_ohm(:,TV1.idx_insp)'));

%% 4) Data selection
% figure; subplot(311); plot(TV.std_RR); subplot(312); plot(TV.std_TV); subplot(313); plot(TV.std_RR.*TV.std_TV);
nan_val = find(isnan(TV1.RR_full));
% last_nan = nan_val(end); %nan value의 마지막 부분
% TV1.RR_isnan =TV1.RR_full(last_nan+1:end); % nan value가 없는 TV.RR 저장

% xlim([-1000 1000] + 473515)
% temp_std_RR_TV = TV1.std_RR.*TV1.std_TV;
% figure; plot(TV1.RR_isnan)
% temp_th_RR = TV1.RR_isnan < 5;
% temp_th_RR = movmean(temp_th_RR, [6000 0]);
% figure; plot(temp_th_RR);

[~,Data.idx_sel] = nanmin(TV1.std_RR.*TV1.std_TV);
Data.idx_sel = Data.idx_sel * op.fs;
% Data.idx_sel = Data.idx_sel - round(60*op.fs/TV1.RR_full(Data.idx_sel));
% if Data.idx_sel < last_nan
% %     Data1.idx_sel = Data1.idx_sel + round(0.5*(60*op.fs/TV.RR_isnan(Data1.idx_sel))); %NAN VALUE 제외하고 시작점 찾기
%     Data.idx_sel = Data.idx_sel+last_nan:Data.idx_sel+last_nan+op.find_ws-1; % nan value로 제외된 숫자 더해서 구간 설정
% else
%     Data1.idx_sel = Data1.idx_sel + round(0.5*(60*op.fs/TV.RR_full(Data1.idx_sel)));
    Data.idx_sel = Data.idx_sel:Data.idx_sel+op.find_ws-1;
% end
% figure; plot(Data.sigma_sum(Data.idx_sel))

if Data.idx_sel(end) > length(TV1.RR_full)
    if Data.idx_sel(end) > length(Data.pixel_sum) % EITdata를 초과할 경우
        Data.mean_RR = nanmean(TV1.RR_full(Data.idx_sel(1):end))/60; % /60: BPM 2 BPS
        Data.std2_RR = 2*nanstd(TV1.RR_full(Data.idx_sel(1):end))/60; % /60: BPM 2 BPS
        TV1.sel_pixel_sum = zeros(op.find_ws,1); % zero value 생성
        EIT_pixel=Data.pixel_sum(Data.idx_sel(1):end); %1~6000중 해당되는 데이터만 추가 
        EIT_pixel = EIT_pixel-mean(EIT_pixel); 
        TV1.sel_pixel_sum(1:length(EIT_pixel)) = EIT_pixel; %zero padding
    else
        Data.mean_RR = nanmean(TV1.RR_full(Data.idx_sel(1):end))/60;
        Data.std2_RR = 2*nanstd(TV1.RR_full(Data.idx_sel(1):end))/60;
        TV1.sel_pixel_sum = Data.pixel_sum(Data.idx_sel);
        TV1.sel_pixel_sum= TV1.sel_pixel_sum-mean(TV1.sel_pixel_sum);
    end
else
    Data.mean_RR = nanmean(TV1.RR_full(Data.idx_sel))/60;
    Data.std2_RR = 2*nanstd(TV1.RR_full(Data.idx_sel))/60;
    TV1.sel_pixel_sum = Data.pixel_sum(Data.idx_sel);
    TV1.sel_pixel_sum= TV1.sel_pixel_sum-mean(TV1.sel_pixel_sum);
end
clear nan_val last_nan;

%% 5) Calc Lung Component
[result.FFT_x_freq,result.FFT_sel] = FxPlotFFT(TV1.sel_pixel_sum,100); xlim([0 5]); ylim([0 1e2]);
result.FFT_sel = abs(result.FFT_sel);
result.FFT_sel = result.FFT_sel./sum(result.FFT_sel); 

if Data.mean_RR-Data.std2_RR<0
    result.sel_start = 1;
    value = [Data.mean_RR Data.std2_RR];
else
   result.sel_start = round((Data.mean_RR-Data.std2_RR)*length(result.FFT_x_freq)/op.fs*2);
end

Lung_1 =(result.sel_start:round((Data.mean_RR+Data.std2_RR)*length(result.FFT_x_freq)/op.fs*2)+1);
Lung_2 =(round((Data.mean_RR*2-Data.std2_RR)*length(result.FFT_x_freq)/op.fs*2):round((Data.mean_RR*2+Data.std2_RR)*length(result.FFT_x_freq)/op.fs*2)+1);
Lung_3 =(round((Data.mean_RR*3-Data.std2_RR)*length(result.FFT_x_freq)/op.fs*2):round((Data.mean_RR*3+Data.std2_RR)*length(result.FFT_x_freq)/op.fs*2)+1);

% lung f0
Lung_co_num = unique(Lung_1);
result.freq_L_f0 = sum(result.FFT_sel(Lung_co_num))*100;

% lung f0 + f1
Lung_co_num = unique([Lung_1 Lung_2]);
result.freq_L_f1 = sum(result.FFT_sel(Lung_co_num))*100;

% lung f0 + f1 + f2
Lung_co_num = unique([Lung_1 Lung_2 Lung_3]);
result.freq_L_f2 = sum(result.FFT_sel(Lung_co_num))*100;

clear Lung_1 Lung_2 Lung_3 Lung_co_num result1.sel_start;

%% 6) drift
op.idx_full = 1:op.ch^2;
op.idx_full(FxEIT_mask(16)) = [];

Data.idx_selv = NaN(256,1);
Data.idx_selv(op.idx_full) = (1:208);
Data.idx_selv = reshape(Data.idx_selv,16,16);
for i = 1:16
    Data.idx_selv(i,:) = circshift(Data.idx_selv(i,:), [0 -i+2]);
end
Data.idx_selv(:,1:3) = [];

ft = fittype( 'poly1' );
for i = 1:16
    Data.EIT_v_ohm_ch_inj(i,:) = mean(Data.EIT_v_ohm((i-1)*13+1:i*13,:));
    temp_1 = TV1.idx_exp - min(TV1.idx_exp);
    temp_2 = Data.EIT_v_ohm_ch_inj(i,TV1.idx_exp);
    [xData, yData] = prepareCurveData( temp_1, temp_2 );
    [fitresult, ~] = fit( xData, yData, ft );
    result.drift_ch_inj(i,1) = fitresult.p1/fitresult.p2*60*op.fs*100; % by 1 min change (%/min)
    result.drift_ch_inj_fitline(i,:) = [min(TV1.idx_exp) max(TV1.idx_exp) fitresult.p1*(min(TV1.idx_exp))+fitresult.p2  fitresult.p1*(max(TV1.idx_exp))+fitresult.p2];
    clear temp_1 temp_2 xData yData fitresult;
    
    Data.EIT_v_ohm_ch_meas(i,:) = mean(Data.EIT_v_ohm(Data.idx_selv(i,:),:));
    temp_1 = TV1.idx_exp - min(TV1.idx_exp);
    temp_2 = Data.EIT_v_ohm_ch_meas(i,TV1.idx_exp);
    [xData, yData] = prepareCurveData( temp_1, temp_2 );
    [fitresult, ~] = fit( xData, yData, ft );
    result.drift_ch_meas(i,1) = fitresult.p1/fitresult.p2*60*op.fs*100; % by 1 min change (%/min)
    result.drift_ch_meas_fitline(i,:) = [min(TV1.idx_exp) max(TV1.idx_exp) fitresult.p1*(min(TV1.idx_exp))+fitresult.p2  fitresult.p1*(max(TV1.idx_exp))+fitresult.p2];
    clear temp_1 temp_2 xData yData fitresult;
end
clear ft
figure;
bar(abs([result.drift_ch_inj result.drift_ch_meas]));

%% 7) Elec Movement
clear indx_near_Erod_single
for i = 1:length(imdl.fwd_model.electrode)
    indx_near_Erod_single{i} = func_get_near_bd_index(imdl.fwd_model.elems, imdl.fwd_model.electrode(i).nodes);
end
% find bnd elem idx
mm = size(imdl.fwd_model.elems,1);
indx_near_bd=func_get_near_bd_index(imdl.fwd_model.elems,imdl.fwd_model.BndIndex);
indx_interior = setdiff(1:mm,indx_near_bd);

% calc S mat
Nskip = 0;
Sensitivity = func_sensitivity_skip(imdl.fwd_model.elems,imdl.fwd_model.nodes,[imdl.fwd_model.electrode.nodes]',ones(mm,1),Nskip,true);
Sensitivity(FxEIT_mask(16),:)=[];

% calc ProjMatMotion
regulMotion = 1e-3;
Proj_MatMotion=func_bd_artifact_elim_new(Sensitivity,indx_interior,regulMotion);
Proj_Mat_ATF = eye(208)-Proj_MatMotion;

bdATFtmp = Proj_Mat_ATF * (Data.EIT_v_ohm(:,Data.idx_sel)-Data.EIT_v_ohm(:,Data.idx_sel(1)));
for i=1:size(indx_near_Erod_single,2)
   bdATF(i,:) = mean(imdl.solve_use_matrix.RM(indx_near_Erod_single{i},:)*bdATFtmp,1);
end
result.sd_bdATF = std(bdATF');
result.bdATF = bdATF;

clear indx_near_Erod_single indx_near_bd indx_interior Sensitivity Nskip mm regulMotion Proj_MatMotion Proj_Mat_ATF bdATF bdATFtmp;

%% disp
figure;
set(gcf,'position',[383,69,757,887]);

% CI
subplot(5,2,2);
bar((result.CI_avg)/1000);
xlabel('Ch','Fontweight','bold');
ylabel('Z_C (k\Omega)','Fontweight','bold');
ylim([0 2]); set(gca,'xtick',1:3:16);
title('Contact impedance','Fontweight','bold');

% RE
subplot(5,2,4);
bar(result.RE_16ch);
xlabel('Ch','Fontweight','bold');
ylabel('RE (%)','Fontweight','bold');
ylim([0 100]); set(gca,'xtick',1:3:16);
title('Reciprocity error','Fontweight','bold');

% Drift I
subplot(5,2,6);
bar(abs(result.drift_ch_inj));
xlabel('Ch_{inj}','Fontweight','bold');
ylabel('Drift (%/min)','Fontweight','bold');
ylim([0 3]); set(gca,'xtick',1:3:16);
title('Drift_{inj}','Fontweight','bold');

% Drift M
subplot(5,2,8);
bar(abs(result.drift_ch_meas));
xlabel('Ch_{meas}','Fontweight','bold');
ylabel('Drift (%/min)','Fontweight','bold');
ylim([0 3]); set(gca,'xtick',1:3:16);
title('Drift_{meas}','Fontweight','bold');

% Movement
subplot(5,2,10);
bar(result.sd_bdATF);
xlabel('Ch','Fontweight','bold');
ylabel('Movement (AU)','Fontweight','bold');
% ylim([0 1]); 
set(gca,'xtick',1:3:16);
title('Electrode Movement','Fontweight','bold');

% Raw V
t = subplot(5,2,[1 3]);
t.Position = [0.13 0.615 0.33 0.30];
plot((1:length(Data.EIT_v(:,Data.idx_sel)))/op.fs/60,Data.EIT_v_ohm(:,Data.idx_sel)','linewidth',2); temp_ylim = ylim; 
xlim([-inf 3000/op.fs/60]); % ylim([0 inf]);
ylim([0.1*temp_ylim(2) 0.25*temp_ylim(2)]);
xlabel('Time (min)','Fontweight','bold'); ylabel('Voltage (V)','Fontweight','bold');
title('Raw data','Fontweight','bold');

% freq
t = subplot(5,2,5);
t.Position = [0.13 0.38 0.33 0.15];
bar(result.FFT_x_freq,result.FFT_sel,'k'); hold on; 
bar(result.FFT_x_freq(result.sel_start:round((Data.mean_RR+Data.std2_RR)*length(result.FFT_x_freq)/op.fs*2)+1),result.FFT_sel(result.sel_start:round((Data.mean_RR+Data.std2_RR)*length(result.FFT_x_freq)/op.fs*2)+1),'r');
bar(result.FFT_x_freq(round((Data.mean_RR*2-Data.std2_RR)*length(result.FFT_x_freq)/op.fs*2):round((Data.mean_RR*2+Data.std2_RR)*length(result.FFT_x_freq)/op.fs*2)+1),result.FFT_sel(round((Data.mean_RR*2-Data.std2_RR)*length(result.FFT_x_freq)/op.fs*2):round((Data.mean_RR*2+Data.std2_RR)*length(result.FFT_x_freq)/op.fs*2)+1),'r');
bar(result.FFT_x_freq(round((Data.mean_RR*3-Data.std2_RR)*length(result.FFT_x_freq)/op.fs*2):round((Data.mean_RR*3+Data.std2_RR)*length(result.FFT_x_freq)/op.fs*2)+1),result.FFT_sel(round((Data.mean_RR*3-Data.std2_RR)*length(result.FFT_x_freq)/op.fs*2):round((Data.mean_RR*3+Data.std2_RR)*length(result.FFT_x_freq)/op.fs*2)+1),'r');
xlim([0 2]); temp_ylim = ylim; ylim([0 temp_ylim(2)*0.5]); 
xlabel('Freq (Hz)','Fontweight','bold'); ylabel('Power (ratio)','Fontweight','bold');
title('Freq Spectrum','Fontweight','bold');

% table
result.RowName = {'L f0','L f0+f1','L f0+f1+f2','RE min','RE avg','RE max','Z insp','Z exp','CI avg','Drift inj','Drift meas'};
result.ColumnName = {'Value','Unit'};
result.table_data = {result.freq_L_f0, result.freq_L_f1, result.freq_L_f2, ...
                     result.RE_min, result.RE_avg, result.RE_max ...
                     result.Z_insp, result.Z_exp, ...
                     nanmean(result.CI_avg) ...
                     round(nanmean(result.drift_ch_inj),4), round(nanmean(result.drift_ch_meas),4)}';
result.table_unit = {'%', '%', '%', '%', '%', '%','ohm','ohm','ohm','%/min','%/min'}';        
t = uitable('Data',[result.table_data, result.table_unit],'RowName',result.RowName,'ColumnName',result.ColumnName);
t.FontSize = 8;
t.Units = 'Normalized';
t.Position = [0.10 0.065 0.361 0.263];


end