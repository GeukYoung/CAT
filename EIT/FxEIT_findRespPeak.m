function [tidal,Detrend_data] = FxEIT_findRespPeak(data,time)
% input
%   data : EIT data 
%   time : EIT time vector
% output
%   tidal : tidal volume & time info

% if nargin < 2
%     time = 1:length(data);
% end
% 
% fs = 1/mean(diff(time)); % EIT framerate
if size(time,1) == 1
    time = time';
end
if size(data,2) == 1
    data = data';
end

fs = 1/mean(diff(time));

%% detrend data  
N = 1; 
% Wn = 0.8;
Wn = 5;
[b,a] = butter(N, 2*Wn/fs,'low');
Detrend_data = filtfilt(b,a,data);
% Detrend_data = data;
dd_data = diff(diff(Detrend_data)); % to find curve shape
%% find insp/exp peak
[~,locs_temp1] = findpeaks(Detrend_data,...
    'MinPeakDistance',round(fs*0.3)); % find vally
[~,locs_temp2] = findpeaks(-Detrend_data,...
    'MinPeakDistance',round(fs*0.3)); % find peak

% figure;
% plot(Detrend_data); hold on;
% plot(locs_temp1,Detrend_data(locs_temp1),'rv','MarkerFaceColor',[1 0 0],'MarkerSize',5);
% plot(locs_temp2,Detrend_data(locs_temp2),'b^','MarkerFaceColor',[0 0 1],'MarkerSize',5);
locs_Peak = [locs_temp1 locs_temp2]';
[locs_Peak,idx_sort] = sort(locs_Peak);

sign_dd_data = [-ones(size(locs_temp1)) ones(size(locs_temp2))]';
sign_dd_data = -sign_dd_data(idx_sort);

% locs_Peak = sort([locs_temp1 locs_temp2]); % combine both peak data
% locs_Peak(locs_Peak>length(dd_data)) = [];
amp_Peak = data(locs_Peak); % amplitude at peak point
% sign_dd_data = sign(dd_data(locs_Peak)); % to find curve shape

% figure;
% plot(sign(dd_data));
% [a, ~] = find(sign_dd_data==1); % find first exp point
% if a(1) > 1 % exp is start
%     locs_Peak(1:a-1) = [];
%     amp_Peak(1:a-1) = [];
%     sign_dd_data(1:a-1) = [];
% end

% figure;
% plot(Detrend_data); hold on;
% plot(locs_Peak,Detrend_data(locs_Peak),'v');

i = 1;
exp_flag = 1;
while(exp_flag)
    if sign_dd_data(1) ~= -1
        sign_dd_data(1) = [];
        locs_Peak(1) = [];
        amp_Peak(1) = [];
        i = i + 1;
    else
        exp_flag = 0;
    end
end
%% insp/exp extraction 
cnt_idx = 2;
sign_curve = 1;
temp_idx_locs = zeros(length(locs_Peak),1);
temp_idx_locs(1) = 1;
for i = 2:length(locs_Peak)-1
    if sign_dd_data(i) ~= sign_curve
        if sign_dd_data(i) == 1
            if amp_Peak(i) > amp_Peak(i-1)
                temp_idx_locs(cnt_idx-1) = i;
            end
        elseif sign_dd_data(i) == -1
            if amp_Peak(i) < amp_Peak(i-1)
                temp_idx_locs(cnt_idx-1) = i;
            end
        end
    else
        temp_idx_locs(cnt_idx) = i;
        cnt_idx = cnt_idx + 1;
        sign_curve = -sign_curve;
    end
end
temp_idx_locs(temp_idx_locs==0) = [];
cnt_idx = cnt_idx-1;

if mod(cnt_idx,2) ~= 0
    cnt_idx = cnt_idx - 1;
end

figure;
plot(Detrend_data); hold on;
plot(locs_Peak(temp_idx_locs(2:2:end)),Detrend_data(locs_Peak(temp_idx_locs(2:2:end))),'rv','MarkerFaceColor',[1 0 0],'MarkerSize',5);
plot(locs_Peak(temp_idx_locs(1:2:end)),Detrend_data(locs_Peak(temp_idx_locs(1:2:end))),'b^','MarkerFaceColor',[0 0 1],'MarkerSize',5);

idx_locs = zeros(cnt_idx,1);
locs_Peak2 = zeros(cnt_idx,1);

idx_locs = temp_idx_locs(1:cnt_idx);
locs_Peak2 = locs_Peak(idx_locs);

tidal_th = 0.3;
idx_rm = [];
for i = 2:(length(locs_Peak2)/2)
    pre_tidal = abs(Detrend_data(locs_Peak2(2*(i-1)-1)) - Detrend_data(locs_Peak2(2*(i-1))));
    cur_tidal = abs(Detrend_data(locs_Peak2(2*i-1)) - Detrend_data(locs_Peak2(2*i)));
    if (tidal_th*pre_tidal) > cur_tidal
        idx_rm = [idx_rm i];
    end
end

if length(idx_rm) > 1
    locs_Peak2([(2*idx_rm-1) (2*idx_rm)]) = [];
end
    
% fig = figure; set(fig,'units','normalized' ,'outerposition' ,[0.15 0.20 0.36 0.31]); 
% for i = 1:(length(locs_Peak2)/2)-1
%     hold on; plot(([locs_Peak2(2*i-1):locs_Peak2(2*i)])/fs ...
%         ,Detrend_data(locs_Peak2(2*i-1):locs_Peak2(2*i)),'k', 'LineWidth',2);
%    hold on; plot(([locs_Peak2(2*i):locs_Peak2(2*i+1)])/fs ...
%         ,Detrend_data(locs_Peak2(2*i):locs_Peak2(2*i+1)),'Color',[0.2 0.9 0.5], 'LineWidth',2);
% end
% hold on; plot((locs_Peak2(1:2:end))/fs,Detrend_data(locs_Peak2(1:2:end)),'b^', ...
%      'MarkerFaceColor',[0 0 1],'MarkerSize',5);
% hold on; plot((locs_Peak2(2:2:end))/fs,Detrend_data(locs_Peak2(2:2:end)),'rv', ...
%      'MarkerFaceColor',[1 0 0],'MarkerSize',5);
% set(gca,'ytick',[]); xlabel('Time(s)','fontsize',8); ylabel('Volume (AU)','fontsize',10);

tidal.idx_exp = zeros(length(locs_Peak2),1);
tidal.idx_insp = zeros(length(locs_Peak2),1);
tidal.t_exp = zeros(length(locs_Peak2),1);
tidal.t_insp = zeros(length(locs_Peak2),1);
tidal.volume_insp = zeros(length(locs_Peak2),1);
tidal.volume_exp = zeros(length(locs_Peak2),1);
tidal.volume_residual = zeros(length(locs_Peak2),1);
tidal.interval = zeros(length(locs_Peak2),1);
tidal.td_insp = zeros(length(locs_Peak2),1);
tidal.td_exp = zeros(length(locs_Peak2),1);
tidal.td_tidal = zeros(length(locs_Peak2),1);
tidal.IE = zeros(length(locs_Peak2),1);
tidal.flow_insp = zeros(length(locs_Peak2),1);
tidal.flow_exp = zeros(length(locs_Peak2),1);
tidal.RR = zeros(length(locs_Peak2),1);

tidal.idx_exp = locs_Peak2(1:2:end)';
tidal.idx_insp = locs_Peak2(2:2:end)';
tidal.t_exp = time(tidal.idx_exp);
tidal.t_insp = time(tidal.idx_insp);
tidal.volume_insp = abs(Detrend_data(tidal.idx_insp) - Detrend_data(tidal.idx_exp))';
tidal.volume_exp = abs(Detrend_data(tidal.idx_insp(1:end-1)) - Detrend_data(tidal.idx_exp(2:end)))';
tidal.volume_residual = tidal.volume_insp(1:end-1) - tidal.volume_exp';
tidal.interval = [0; (time(tidal.idx_insp(2:end))-time(tidal.idx_insp(1:end-1)))];
tidal.td_insp = time(tidal.idx_insp) - time(tidal.idx_exp);
tidal.td_exp = time(tidal.idx_exp(2:end)) - time(tidal.idx_insp(1:end-1));
tidal.td_tidal = tidal.td_insp(1:end-1) + tidal.td_exp;
tidal.IE = tidal.td_insp(1:end-1)./tidal.td_tidal;
tidal.flow_insp = tidal.volume_insp./tidal.td_insp;
tidal.flow_exp = tidal.volume_exp./tidal.td_exp;
tidal.RR = 60./(time(tidal.idx_insp(2:end))-time(tidal.idx_insp(1:end-1)));
end

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
% % 	hold on; plot(([locs_Peak(2*i):locs_Peak(2*i+1)])/fs ...
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