function [Score] = FxPSG_Scoring(data,time,PSG,tidal,time_section)
% input
%   data : EIT data vector
%   time : EIT time vector
%   PSG : PSG data
%   tidal : tidal volume & time info
%   time_section : movement segmentation
% output
%   Score : result of scoring

if nargin < 5
    time_section = time(tidal.idx_insp(end));
else
    time_section = [time_section time(tidal.idx_insp(end))];
end

% Event rules
rule.Apnea.time = 10; % 10s duration
rule.Hypop.amp_th = 0.7; % 30% signal Drop -> Hypopnea
rule.Apnea.amp_th = 0.20; % 80% signal Drop -> Apnea
rule.Apnea.cen_th = 0.02; % 98% signal Drop -> Central Apnea
rule.Hypop.desat_time = 20; % time extand 20s length
rule.Hypop.amp_snor = 0.00002; % snoring sound threshold (20u)
rule.RERA.time = 10; % flattening or effort 10s
rule.Epoch = 10*30;

% Event Definition
event.Normal = 0;
event.Obstructive_Apnea = 1;
event.Mixed_Apnea = 2;
event.Central_Apnea = 3;
event.Obstructive_Hypopnea = 4;
event.Central_Hypopnea = 5;
event.RERA = 6;
event.tag = {'Normal','ObstructiveApnea','MixedApnea','CentralApnea', ...
    'ObstructiveHypopnea','CentralHypopnea','RERA'};

% parameter initial
state = 0;
sEvent = 0;
cnt_Event = 1;
cnt_Section = 1;
time_apnea = 0;
time_hypop = 0;
time_central = 0;
amp_normalTV = 0; 
idx_event = zeros(1,length(tidal.volume_insp));
idx_ref = 1;
tidal.time_insp = time(tidal.idx_insp);

for cntTidal = 2:length(tidal.volume_insp)
    % Scoring segmentation for movement
    if tidal.time_insp(cntTidal) > time_section(cnt_Section)
        cnt_Section = cnt_Section + 1;
        state = 0; % state reset
        amp_normalTV = 0; 
    end

    switch state
       % Normal
        case 0
           % Check Tidal Amp Drop
            if tidal.volume_insp(cntTidal) < amp_normalTV*rule.Hypop.amp_th
                state = 1; % Event
                sEvent = cntTidal;
           % Check Central Apnea
            elseif time(tidal.idx_exp(cntTidal))-time(tidal.idx_insp(cntTidal-1)) > rule.Apnea.time
                % central
                event_type = event.Central_Apnea;  
                sEvent = cntTidal-1;
                state = 3;
           % NormalTidal Updata
            else
                if tidal.time_insp(cntTidal) < time_section(cnt_Section) - rule.Epoch
                    tidal_sel = tidal.time_insp;
                    tidal_sel(tidal_sel < tidal.time_insp(cntTidal)) = 0;
                    tidal_sel(tidal_sel > tidal.time_insp(cntTidal) + rule.Epoch) = 0;
                    tidal_sel(tidal_sel > 0) = 1;
                    tidal_sel = tidal_sel .* (1:length(tidal_sel))';
                    tidal_sel(tidal_sel == 0) = [];
                    
                    temp = tidal.volume_insp(tidal_sel);
                    temp(:,2) = tidal_sel;
                    temp = sortrows(temp);
                    ref_th = round(length(temp)*0.1);
                    [amp_normalTV, idx_Tidal_normal] = FxOptimal_Lcurve(temp(1+ref_th:end-ref_th,:));
                    idx_ref = idx_Tidal_normal;
                end
            end
            
       % Event
        case 1
           % Check Amp Recorver
            if tidal.volume_insp(cntTidal) < amp_normalTV*rule.Hypop.amp_th
                % need adding watchdog
                % state = 0;
            else
               % Check During Time
                if time(tidal.idx_exp(cntTidal))-time(tidal.idx_exp(sEvent)) > rule.Apnea.time
                 % Apnea Classification
                   % Apnea
                    event_volumes = tidal.volume_insp(sEvent:cntTidal);
                    event_intervals = tidal.interval(sEvent:cntTidal);
                    time_central = sum(event_intervals(event_volumes<amp_normalTV*rule.Apnea.cen_th));
                    time_apnea = sum(event_intervals(event_volumes<amp_normalTV*rule.Apnea.amp_th));
                    time_hypop = sum(event_intervals(event_volumes<amp_normalTV*rule.Hypop.amp_th)) - time_apnea;
                                        
                    if time_apnea > rule.Apnea.time
                        if time_central > rule.Apnea.time
                            % central
                            event_type = event.Central_Apnea;  
                            state = 3;
                        elseif time_central > rule.Apnea.time/2
                            % Mixed
                            event_type = event.Mixed_Apnea;  
                            state = 3;
                        else
                            % Apnea
                            event_type = event.Obstructive_Apnea;  
                            state = 3;
                        end
                   % Hypopnea
                    else
                       % Check SpO2 DeSat                             
                        temp_SpO2 = sum(PSG{7}.Desat_flag(FxFindIdx(PSG{7}.timeinfo,time(tidal.idx_exp(sEvent))): ...
                             FxFindIdx(PSG{7}.timeinfo,time(tidal.idx_exp(cntTidal))+rule.Hypop.desat_time)));
                        temp_SpO2(temp_SpO2>0) = 1;
                        if temp_SpO2 == 1
                            temp_Snor = PSG{9}.data(FxFindIdx(PSG{9}.timeinfo,time(tidal.idx_exp(sEvent))) : ...
                                FxFindIdx(PSG{9}.timeinfo,time(tidal.idx_exp(cntTidal))));
                            if max(temp_Snor) > rule.Hypop.amp_snor % Snoring
                                event_type = event.Obstructive_Hypopnea;  
                                state = 3;
                            else % can't detect effort yet
                                event_type = event.Obstructive_Hypopnea;  
                                state = 3;
%                           elseif 2 % flattening
%                               event_type = event.Obstructive_Hypopnea;  
%                               state = 3;
%                           elseif 3 % effort
%                               event_type = event.Obstructive_Hypopnea;  
%                               state = 3;
%                           else
%                               event_type = event.Central_Hypopnea;  
%                               state = 3;
                            end
                        else
                            event_type = event.Normal;  
                            state = 3;
                        end
                    end
                else
                    event_type = event.Normal;  
                    state = 3;
                end
            end
    end
    
   % Save result
    if state == 3
        if event_type ~= event.Normal
            Score.event{cnt_Event}.time_start = time(tidal.idx_exp(sEvent));
            Score.event{cnt_Event}.time_end = time(tidal.idx_exp(cntTidal));
            Score.event{cnt_Event}.duration = Score.event{cnt_Event}.time_end - Score.event{cnt_Event}.time_start;
            Score.event{cnt_Event}.idx_type = event_type;
            Score.event{cnt_Event}.idx_tidal = sEvent:cntTidal;
            Score.event{cnt_Event}.idx_ref = idx_ref;
            Score.event{cnt_Event}.lv_th = amp_normalTV;
            Score.event{cnt_Event}.type = event.tag{event_type};
            Score.event{cnt_Event}.time_apnea = time_apnea;
            Score.event{cnt_Event}.time_hypop = time_hypop;
            Score.event{cnt_Event}.time_central = time_central;
            
            if event_type == event.Central_Apnea
                Score.event{cnt_Event}.time_start = time(tidal.idx_insp(cntTidal));
                Score.event{cnt_Event}.time_end = time(tidal.idx_exp(cntTidal));
                Score.event{cnt_Event}.duration = Score.event{cnt_Event}.time_end - Score.event{cnt_Event}.time_start;
                Score.event{cnt_Event}.time_apnea = 0;
                Score.event{cnt_Event}.time_hypop = 0;
                Score.event{cnt_Event}.time_central = Score.event{cnt_Event}.duration;
            end
            idx_event(sEvent:cntTidal) = event_type.*ones(1,cntTidal-sEvent+1);
            cnt_Event = cnt_Event + 1;
        end
        % init
        state = 0;
        eEvent = cntTidal;
%         temp = tidal.volume_insp(cntTidal:cntTidal+rule.Tidal.windowsize);
%         Tidal_normal = max([temp(temp<max(temp));Tidal_normal]); clear temp;
    end
    
    disp(strcat(num2str(cntTidal),'/',num2str(length(tidal.volume_insp))));
end

% Makeup score result
for i = 1:length(Score.event)
    temp_event(i) = Score.event{i}.idx_type;
end

Score.Obstructive_Apnea = sum(temp_event==event.Obstructive_Apnea);
Score.Mixed_Apnea = sum(temp_event==event.Mixed_Apnea);
Score.Central_Apnea = sum(temp_event==event.Central_Apnea);
Score.Obstructive_Hypopnea = sum(temp_event==event.Obstructive_Hypopnea);
Score.Central_Hypopnea = sum(temp_event==event.Central_Hypopnea);

disp(['Obstructive_Apnea : ' ,num2str(Score.Obstructive_Apnea)]);
disp(['Mixed_Apnea : ' ,num2str(Score.Mixed_Apnea)]);
disp(['Central_Apnea : ' ,num2str(Score.Central_Apnea)]);
disp(['Obstructive_Hypopnea : ' ,num2str(Score.Obstructive_Hypopnea)]);
disp(['Central_Hypopnea : ' ,num2str(Score.Central_Hypopnea)]);
AHI = length(Score.event)/((time(end)-time(1))/60^2);
disp(['AHI : ' ,num2str(AHI)]);

idxPSG.EKG = 1;
idxPSG.NasalPressure = 2;
idxPSG.Flow_CU = 3;
idxPSG.Thermistor = 4;
idxPSG.Thorax = 5;
idxPSG.Abdomen = 6;
idxPSG.SpO2 = 7;
idxPSG.Position_PDS = 8;
idxPSG.Activity_CU = 9;
idxPSG.A1 = 10;
idxPSG.A2 = 11;
idxPSG.C3 = 12;
idxPSG.C4 = 13;
idxPSG.F3 = 14;
idxPSG.F4 = 15;
idxPSG.O1 = 16;
idxPSG.O2 = 17;
idxPSG.Snoring_Sensor = 18;
idxPSG.Snore_CU = 19;
%% Display result
fs = 50;
% fig = figure; set(fig,'units','normalized' ,'outerposition' ,[0.15 0.20 0.36 0.31]); 

clear fig;
figure(); set(gcf,'units','normalized','outerposition',[0 0 0.7 0.82]);
n_lowsampling = 1;
n_signal = 7;

% EIT
fig{1} = subplottight(n_signal,1,1,2);
plot(time(1:n_lowsampling:end), ...
    data(1:n_lowsampling:end),'LineWidth',2);
set(fig{1},'ytick',[]); set(fig{1},'xtick',[]);
ylabel('Tidal Volume','fontsize',15,'Rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','right');

max_tidal = max(data) * 1.1;
min_tidal = min(data);
% max_tidal = 2*10^6;
% min_tidal = 0;
for i = 1:length(Score.event)
    hold on;
    vert = [Score.event{i}.time_start min_tidal; ...
        Score.event{i}.time_start max_tidal; ...
        Score.event{i}.time_end max_tidal; ...
        Score.event{i}.time_end min_tidal];
    fac = [1 2 3 4];
    patch('Faces',fac,'Vertices',vert,'FaceColor','flat','FaceVertexCData',Score.event{i}.idx_type,'FaceAlpha',.9); colormap jet; caxis([0 5]); 
%     plot([Score.event{i}.time_start Score.event{i}.time_start],[data(tidal.idx_exp(Score.event{i}.idx_tidal(1))) data(tidal.idx_exp(Score.event{i}.idx_tidal(1)))+Score.event{i}.lv_th],'linewidth',5,'color','w');
%     text(Score.event{i}.time_start,min_tidal,['#',num2str(i),' (',num2str(round(Score.event{i}.duration)),'s)'],'fontsize',12);
end

fig{2} = subplottight(n_signal,1,2,2);
plot(PSG{idxPSG.NasalPressure}.timeinfo(1:n_lowsampling:end), ...
    PSG{idxPSG.NasalPressure}.filt_data(1:n_lowsampling:end),'LineWidth',2,'Color',[0.85 0.33 0.1]);
set(fig{2},'ytick',[]); set(fig{2},'xtick',[]);
ylabel('Nasal Pressure','fontsize',15,'Rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','right');


fig{3} = subplottight(n_signal,1,3,2);
plot(PSG{idxPSG.Thermistor}.timeinfo(1:n_lowsampling:end), ...
    PSG{idxPSG.Thermistor}.filt_data(1:n_lowsampling:end),'LineWidth',2,'Color',[0.87 0.49 0]);
set(fig{3},'ytick',[]); set(fig{3},'xtick',[]);
ylabel('Thermistor','fontsize',15,'Rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','right');

fig{4} = subplottight(n_signal,1,4,2); hold on;
plot(PSG{idxPSG.Thorax}.timeinfo(1:n_lowsampling:end), ...
    3*PSG{idxPSG.Thorax}.filt_data(1:n_lowsampling:end)+2e-4,'LineWidth',2,'Color',[0 0.5 0]);
plot(PSG{idxPSG.Abdomen}.timeinfo(1:n_lowsampling:end), ...
    PSG{idxPSG.Abdomen}.filt_data(1:n_lowsampling:end)-2e-4,'LineWidth',2,'Color',[0.75 0 0.75]);
set(fig{4},'ytick',[]); set(fig{4},'xtick',[]);
ylabel({'Thorax';'Abdomen'},'fontsize',15,'Rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','right');

fig{5} = subplottight(n_signal,1,5,2);
hold on;
for cnt = 1:length(PSG{idxPSG.SpO2}.DesatEvent)
    vert = [PSG{idxPSG.SpO2}.timeinfo(PSG{idxPSG.SpO2}.DesatEvent(cnt,1)) 80; ...
        PSG{idxPSG.SpO2}.timeinfo(PSG{idxPSG.SpO2}.DesatEvent(cnt,1)) 100; ...
        PSG{idxPSG.SpO2}.timeinfo(PSG{idxPSG.SpO2}.DesatEvent(cnt,2)) 100; ...
        PSG{idxPSG.SpO2}.timeinfo(PSG{idxPSG.SpO2}.DesatEvent(cnt,2)) 80];
    fac = [1 2 3 4];
    patch('Faces',fac,'Vertices',vert,'FaceColor','red','FaceAlpha',.7);
    %         colormap jet; caxis([-1.5 1.5]);
end
plot(PSG{idxPSG.SpO2}.timeinfo(1:n_lowsampling:end), ...
    PSG{idxPSG.SpO2}.filt_data(1:n_lowsampling:end),'LineWidth',2);
set(fig{5},'Tickdir','in');
set(fig{5},'xtick',[]);
ylim([83 100]);
ylabel('SpO2','fontsize',15,'Rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','right');


fig{6} = subplottight(n_signal,1,6,2);
plot(PSG{idxPSG.Snore_CU}.timeinfo(1:n_lowsampling:end), ...
    PSG{idxPSG.Snore_CU}.data(1:n_lowsampling:end),'k','LineWidth',2);
set(fig{6},'ytick',[]); set(fig{6},'xtick',[time(1):10:time(end)]); set(fig{6},'xticklabels',[]);
ylabel('Snore Sound','fontsize',15,'Rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','right');
xlabel('Time(s)','fontsize',14);
% set(fig{9},'xtick',EIT.timeinfo(1):10:EIT.timeinfo(end));
% grid on;

linkaxes([fig{:}],'x');
set(gcf, 'color', 'white');
for i = 1:6
    set(fig{i},'Box','off');
end

% ticklabelinside(fig{5},'y',1);

% subplottight(4,1,[1 2],2);
% max_tidal = max(data);
% min_tidal = min(data);
max_tidal = 2*10^6;
min_tidal = 0;
for i = 1:length(Score.event)
    hold on;
    vert = [Score.event{i}.time_start min_tidal; ...
        Score.event{i}.time_start max_tidal; ...
        Score.event{i}.time_end max_tidal; ...
        Score.event{i}.time_end min_tidal];
    fac = [1 2 3 4];
    patch('Faces',fac,'Vertices',vert,'FaceColor','flat','FaceVertexCData',Score.event{i}.idx_type); colormap jet; caxis([0 5]); 
%     plot([Score.event{i}.time_start Score.event{i}.time_start],[data(tidal.idx_exp(Score.event{i}.idx_tidal(1))) data(tidal.idx_exp(Score.event{i}.idx_tidal(1)))+Score.event{i}.lv_th],'linewidth',5,'color','w');
%     text(Score.event{i}.time_start,min_tidal,['#',num2str(i),' (',num2str(round(Score.event{i}.duration)),'s)'],'fontsize',12);
end

for i = 1:length(tidal.idx_exp)-2
    hold on; plot(time([tidal.idx_exp(i):tidal.idx_insp(i)]) ...
        ,data([tidal.idx_exp(i):tidal.idx_insp(i)]),'k', 'LineWidth',2);
	hold on; plot(time([tidal.idx_insp(i):tidal.idx_exp(i+1)]) ...
        ,data(tidal.idx_insp(i):tidal.idx_exp(i+1)),'Color',[0.2 0.9 0.5], 'LineWidth',2);
%     text(time(tidal.idx_insp(i)),data(tidal.idx_insp(i)), ...
%         [' #',num2str(i),' (',num2str(round((tidal.volume_insp(i)/refVolume(i))*100)/100),')'], ...
%         'fontsize',8);
end
plot(time(tidal.idx_insp),tidal.volume_insp,'-o','LineWidth',1.5); 

hold on; plot(time(tidal.idx_exp),data(tidal.idx_exp),'b^', ...
     'MarkerFaceColor',[0 0 1],'MarkerSize',5);
hold on; plot(time(tidal.idx_insp),data(tidal.idx_insp),'rv', ...
     'MarkerFaceColor',[1 0 0],'MarkerSize',5);
set(gca,'ytick',[]); xlabel('Time(s)','fontsize',8); ylabel('Volume (AU)','fontsize',10);
set(gcf, 'color', 'white');
set(gca, 'Xtick', tidal.idx_insp(1):10:tidal.idx_exp(end));
% figure;
plot(PSG{7}.timeinfo,PSG{7}.Desat_flag*min(data),'r');
plot(PSG{7}.timeinfo,FxFilt_Norm(PSG{7}.data)*min(data),'k');
grid on;
% axis([52234.3639481902,52400.8811288558,8969917.06588464,12530070.3409129]);
hold on; plot(PSG{7}.timeinfo,FxFilt_Norm(PSG{7}.data)); 
plot(time,FxFilt_Norm(data)+1,'r');


%%
% figure;
% subplottight(4,1,[3],2);
% hold off;
% plot(PSG{7}.timeinfo,PSG{7}.filt_data,'LineWidth',1);
% ylabel('Nasal Flow'); ylim([50 100]);
% grid on;
% 
% subplottight(4,1,[3],2);
% plot(time(tidal.idx_insp([1 end])),[0 0],'-k','LineWidth',1.5);
% hold on; plot(time(tidal.idx_insp(1:end-1)),tidal.volume_residual,'-o', 'LineWidth',2);
% ylabel('Residual change');
% grid on;
% 
% subplottight(4,1,[3],2);
% hold off;
% plot(PSG{3}.timeinfo,PSG{3}.norm_data,'LineWidth',1);
% ylabel('Nasal Flow');
% grid on;
% 
% 
% subplottight(4,1,[4],2);
% plot(PSG{5}.timeinfo,PSG{5}.norm_data,'LineWidth',1);
% hold on; plot(PSG{6}.timeinfo,PSG{6}.norm_data,'LineWidth',1);
% ylabel('Chest & Abdomen'); legend({'Chest','Abdomen'});
% grid on;
% 
% % tidal.volume_residual = tidal.volume_insp - tidal.volume_exp
% subplottight(4,1,[1 2],2);
% temp = axis;
% subplottight(4,1,[3],2);
% xlim([temp(1) temp(2)]);
% subplottight(4,1,[4],2);
% xlim([temp(1) temp(2)]);
% 
% subplottight(4,1,[1 2],2);
% temp = axis;
% subplottight(4,1,[3]);
% xlim([temp(1) temp(2)]);
% 
% subplottight(4,1,[4],2);
% 
% [temp,~] = ginput(2);
% disp(['Time diff : ',num2str(diff(temp)),'s']);
% 
% 
% fig = figure; set(fig,'units','normalized' ,'outerposition' ,[0.15 0.20 0.36 0.31]); 
% for i = 1:(length(locs_Peak)/2)-1
%     hold on; plot(([locs_Peak(2*i-1):locs_Peak(2*i)])/fs ...
%         ,data(locs_Peak(2*i-1):locs_Peak(2*i)),'k', 'LineWidth',2);
% 	hold on; plot(([locs_Peak(2*i):locs_Peak(2*i+1)])/fs ...
%         ,data(locs_Peak(2*i):locs_Peak(2*i+1)),'Color',[0.2 0.9 0.5], 'LineWidth',2);
% end
% hold on; plot((locs_Peak(1:2:end))/fs,data(locs_Peak(1:2:end)),'b^', ...
%      'MarkerFaceColor',[0 0 1],'MarkerSize',5);
% hold on; plot((locs_Peak(2:2:end))/fs,data(locs_Peak(2:2:end)),'rv', ...
%      'MarkerFaceColor',[1 0 0],'MarkerSize',5);
% set(gca,'ytick',[]); xlabel('Time(s)','fontsize',8); ylabel('Volume (AU)','fontsize',10);
% xlim([7501 9000]);
