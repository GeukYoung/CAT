function [Beat,idx_En,locs_Rwave] = FxSCG_FPs(SCG,ECG,fs)
ECG(min([length(SCG), length(ECG)]):end) = [];

%% Data Setup
Seismo.raw = SCG; % Acc Z is seismo

rule.seg_R = 0.200; % 200 ms segmentation from R peak
rule.idx_R = round(rule.seg_R * fs);
rule.min_S1rd = round(0.010 * fs); % 10 ms S1 Rdelay threshold (min)
rule.max_S1rd = round(0.160 * fs); % 160 ms S1 Rdelay threshold (max)
rule.min_S2rd = round(0.300 * fs); % 300 ms S2 Rdelay threshold (min)
rule.max_S2rd = round(0.480 * fs); % 480 ms S2 Rdelay threshold (max)
rule.idx_thMC = round(0.05 * fs); % 50 ms minimum distance from ICP
rule.th_AO = 0.5; % AO amp threshold greater than 0.5 * ICPd (0.7 happen missed case)
rule.th_AC = 1.2; % AC amp threshold (D1 > th_AC*(IRP-AC))
rule.th_ICP = 0.9;
rule.idx_thS2si = round(0.1 * fs); % S2 range 60 ms -> 100 ms

%% Part 1. Post process (bp filt & envelope extraction)
n = 1; fc = [5 40];
Seismo.filt = FxEIT_Filter(Seismo.raw,fs,n,fc,1,3);
n = fs/10; b = [1:n n-1:-1:1]/n^2*1.75;
Seismo.envelope = filtfilt(b,1,abs(Seismo.filt)); 

% subplottight(5,1,1);
% plot(Seismo.raw); set(gca,'XTick',[],'YTick',[]); xlim([1.5 1.8]*10^4);
% subplottight(5,1,2);
% plot(Seismo.filt); set(gca,'XTick',[],'YTick',[]); xlim([1.5 1.8]*10^4);
% subplottight(5,1,3); 
% plot(Seismo.filt); hold on;
% plot(Seismo.envelope,'k'); xlim([1.5 1.8]*10^4);
% hold on;
% subplottight(5,1,4); 
% plot(ECG,'r'); xlim([1.5 1.8]*10^4);

%% Part 2. Beat to Beat segmentation using ECG
[locs_Rwave,~,~,~] = FxEIT_findRpeak(ECG,fs);
locs_Rwave = locs_Rwave - 9;
locs_Rwave(1) = [];
temp_seg = locs_Rwave - rule.seg_R*fs + 1;
nBeat = length(temp_seg)-1;
for cntBeat = 1:nBeat
    Beat(cntBeat).idx_raw = [temp_seg(cntBeat):temp_seg(cntBeat+1)-1];
    Beat(cntBeat).ecg = ECG(Beat(cntBeat).idx_raw);
    Beat(cntBeat).seismo = Seismo.filt(Beat(cntBeat).idx_raw);
    Beat(cntBeat).seismo_envelope = Seismo.envelope(Beat(cntBeat).idx_raw);
    Beat(cntBeat).idx_R = rule.idx_R ;
    Beat(cntBeat).flag_n = 0;
    
    [~, locs_Q] = findpeaks(-Beat(cntBeat).ecg);
%     figure; plot(temp); hold on; plot(locs_Qb,-a,'rv');
    locs_Q = locs_Q(locs_Q<rule.idx_R);
    Beat(cntBeat).idx_Q = locs_Q(end);
end
disp(['1) ECG segment num : ',num2str(nBeat)]);

% subplottight(5,1,4); 
% plot(ECG); hold on;
% plot(locs_Rwave,ECG(locs_Rwave),'rv'); xlim([1.5 1.8]*10^4);
% plot(temp_seg,zeros(length(temp_seg),1),'ro');

%% Part 3. S1 & S2 detection
for cntBeat = 1:nBeat
    if Beat(cntBeat).flag_n == 0
        [Amp_Swave, locs_Swave] = findpeaks(Beat(cntBeat).seismo_envelope(rule.seg_R*fs+1:end));
        locs_Swave = locs_Swave + rule.seg_R*fs;
        if length(locs_Swave) < 2 % Exception1. number of peak
            Beat(cntBeat).flag_n = 1;
        elseif Amp_Swave(2) > Amp_Swave(1) % Exception2. S1 > S2
            Beat(cntBeat).flag_n = 1;
        elseif ((locs_Swave(1)-rule.seg_R*fs) < rule.min_S1rd) || ...
                ((locs_Swave(1)-rule.seg_R*fs) > rule.max_S1rd) % Exception3. check S1rd
            Beat(cntBeat).flag_n = 1;
        elseif ((locs_Swave(2)-rule.seg_R*fs) < rule.min_S2rd) || ...
                ((locs_Swave(2)-rule.seg_R*fs) > rule.max_S2rd) % Exception4. check S2rd
            Beat(cntBeat).flag_n = 1;
        else % Stack S1 S2 time
            Beat(cntBeat).idx_S1 = locs_Swave(1);
            Beat(cntBeat).idx_S2 = locs_Swave(2);
        end
    end
end
disp(['2) S1,S2 detection : ',num2str(sum([Beat.flag_n])),'/',num2str(nBeat)]);
% figure;
% plot([Beat.seismo_envelope]); hold on;
% plot([Beat.seismo]);
% for cntBeat = 1:nBeat
%     plot(Beat(cntBeat).idx_S1+Beat(cntBeat).idx(1)-1,Seismo.envelope(Beat(cntBeat).idx_S1+Beat(cntBeat).idx(1)-1),'rv')
%     plot(Beat(cntBeat).idx_S2+Beat(cntBeat).idx(1)-1,Seismo.envelope(Beat(cntBeat).idx_S2+Beat(cntBeat).idx(1)-1),'bv')
% end
% 
% % plot result
% figure;
% subplottight(5,1,1); 
% plot(ECG); hold on;
% plot(locs_Rwave,ECG(locs_Rwave),'rv'); xlim([1.5 1.8]*10^4);
% plot(temp_seg,zeros(length(temp_seg),1),'ro');
% subplottight(5,1,2);
% plot(Seismo.raw); set(gca,'XTick',[],'YTick',[]); xlim([1.5 1.8]*10^4);
% subplottight(5,1,3);
% plot(Seismo.filt); set(gca,'XTick',[],'YTick',[]); xlim([1.5 1.8]*10^4);
% subplottight(5,1,4); 
% plot(abs(Seismo.filt)); hold on;
% plot(Seismo.envelope,'k'); xlim([1.5 1.8]*10^4);
% for cntBeat = 1:nBeat
%     plot(Beat(cntBeat).idx_S1+Beat(cntBeat).idx(1)-1,Seismo.envelope(Beat(cntBeat).idx_S1+Beat(cntBeat).idx(1)-1),'rv')
%     plot(Beat(cntBeat).idx_S2+Beat(cntBeat).idx(1)-1,Seismo.envelope(Beat(cntBeat).idx_S2+Beat(cntBeat).idx(1)-1),'bv')
% end
% subplottight(5,1,5);
% plot(Seismo.filt); set(gca,'XTick',[],'YTick',[]); xlim([1.5 1.8]*10^4); hold on;
% for cntBeat = 1:nBeat
%     plot(Beat(cntBeat).idx_S1+Beat(cntBeat).idx(1)-1,Seismo.envelope(Beat(cntBeat).idx_S1+Beat(cntBeat).idx(1)-1),'rv')
%     plot(Beat(cntBeat).idx_S2+Beat(cntBeat).idx(1)-1,Seismo.envelope(Beat(cntBeat).idx_S2+Beat(cntBeat).idx(1)-1),'bv')
% end

%% Part 4. FPs from S1 signal (ICP MC AO)
cntBeat = 158
cntBeat = 4
for cntBeat = 1:nBeat
    if Beat(cntBeat).flag_n == 0
        S1 = Beat(cntBeat).seismo(rule.idx_R:rule.idx_R + 2*(Beat(cntBeat).idx_S1 - rule.idx_R));
        %% ICP detection
        [Amp_FP, locs_FP] = findpeaks(-S1);
        locs_FP = locs_FP(locs_FP < Beat(cntBeat).idx_S1 - rule.idx_R + 1); % choose only peaks before S1
        %     plot(S1); hold on; plot(locs_FP,S1(locs_FP),'rv');
        
        for cntFP = 1:length(locs_FP)
            if (S1(locs_FP(cntFP)) < min(S1(locs_FP)) * rule.th_ICP)
                Beat(cntBeat).idx_ICP = locs_FP(cntFP) + rule.idx_R - 1; % phase compantation (cuz extract data set)
                Beat(cntBeat).amp_ICP = Amp_FP(cntFP);
            end
        end
        if isempty(Beat(cntBeat).idx_ICP) == true % if do not find ICP => fail flag
            Beat(cntBeat).flag_n = 1; 
        end
        
         %% MC & AO detection
        [~, locs_FP] = findpeaks(S1);
        %     plot(S1); hold on; plot(locs_FP,S1(locs_FP),'bv');
%         figure;
%         plot(Beat(cntBeat).ecg);
%         figure;
%         plot(Beat(cntBeat).seismo);
%         figure;
%         plot(Beat(cntBeat).seismo_envelope);
        
        % MC detection
        locs_MC = locs_FP(locs_FP < Beat(cntBeat).idx_ICP - rule.idx_R + 1); % choose only peaks befor ICP
        if  (locs_MC(end) + rule.idx_R - 1 > Beat(cntBeat).idx_ICP - rule.idx_thMC) ... % check time threshold
                && (S1(locs_MC(end)) > 0) % check amplitude threshold
            Beat(cntBeat).idx_MC = locs_MC(end) + rule.idx_R - 1; % phase compantation (cuz extract data set)
        else
            Beat(cntBeat).flag_n = 1;
        end
        
        % AO detection
        Beat(cntBeat).idx_AO = [];
        locs_AO = locs_FP(locs_FP > Beat(cntBeat).idx_ICP - rule.idx_R + 1); % choose only peaks after ICP
        if isempty(locs_AO) == false
            for cntFP = 1:length(locs_AO)
                if (S1(locs_AO(cntFP)) > Beat(cntBeat).amp_ICP * rule.th_AO) && isempty(Beat(cntBeat).idx_AO)
                    Beat(cntBeat).idx_AO = locs_AO(cntFP) + rule.idx_R - 1; % phase compantation (cuz extract data set)
                end
            end
        end
        if isempty(Beat(cntBeat).idx_AO) == true % if do not find AO => fail flag
            Beat(cntBeat).flag_n = 1;
        end
    end
    disp([num2str(cntBeat)]);
end
disp(['3) ICP,MC,AO detection : ',num2str(sum([Beat.flag_n])),'/',num2str(nBeat)]);

%% Part 5. FPs from S2 signal (IRP AC MO)
for cntBeat = 1:nBeat
    if Beat(cntBeat).flag_n == 0
        try
            S2 = Beat(cntBeat).seismo((Beat(cntBeat).idx_S2-rule.idx_thS2si):(Beat(cntBeat).idx_S2+rule.idx_thS2si));
        catch
            S2 = Beat(cntBeat).seismo((Beat(cntBeat).idx_S2-rule.idx_thS2si):end);
        end
        plot(Beat(cntBeat).seismo);
%         plot(S2);
%         plot(Beat(cntBeat).seismo); hold on;
%         plot(Beat(cntBeat).seismo_envelope);
%         
%         plot(Beat(cntBeat-2).seismo); hold on;
%         plot(Beat(cntBeat-2).seismo_envelope);
        
        %% IRP detection
        [IRPd, IRP] = max(S2);
        Beat(cntBeat).idx_IRP = IRP + Beat(cntBeat).idx_S2 - rule.idx_thS2si - 1; % phase compantation (cuz extract data set)
        Beat(cntBeat).amp_IRPd = IRPd;

        %% AC detection
        [~ , locs_FP] = findpeaks(S2);
        %     figure; plot(S2); hold on; plot(locs_FP,S2(locs_FP),'bv');
        locs_AC = locs_FP(locs_FP < Beat(cntBeat).idx_IRP - Beat(cntBeat).idx_S2 + rule.idx_thS2si);
        for i = 1:length(locs_AC)
            tempD1(i) = S2((Beat(cntBeat).idx_IRP - Beat(cntBeat).idx_S2 + rule.idx_thS2si)) - ...
                min(S2(locs_AC(i):(Beat(cntBeat).idx_IRP - Beat(cntBeat).idx_S2 + rule.idx_thS2si)));
            tempIRP_AC(i) = S2((Beat(cntBeat).idx_IRP - Beat(cntBeat).idx_S2 + rule.idx_thS2si)) - ...
                S2(locs_AC(i));
        end
        
        Beat(cntBeat).idx_AC = [];
        for i = length(locs_AC):-1:1
            if (tempD1(i) > tempIRP_AC(i) * rule.th_AC) && isempty(Beat(cntBeat).idx_AC)
                Beat(cntBeat).idx_AC = locs_AC(i) + Beat(cntBeat).idx_S2 - rule.idx_thS2si - 1; % phase compantation (cuz extract data set)
                Beat(cntBeat).D1 = tempD1(i);
            end
        end
        if isempty(Beat(cntBeat).idx_AC) == true
            Beat(cntBeat).flag_n = 1;
        end
        
        %% MO detection
        [~ , locs_FP] = findpeaks(-S2);
        %     plot(S2); hold on; plot(locs_FP,S2(locs_FP),'bv');
        locs_MO = locs_FP(locs_FP > Beat(cntBeat).idx_IRP - Beat(cntBeat).idx_S2 + rule.idx_thS2si);
        if isempty(locs_MO) == false
            Beat(cntBeat).idx_MO = locs_MO(1) + Beat(cntBeat).idx_S2 - rule.idx_thS2si - 1; % phase compantation (cuz extract data set)
            Beat(cntBeat).D2 = Beat(cntBeat).seismo(Beat(cntBeat).idx_IRP) - Beat(cntBeat).seismo(Beat(cntBeat).idx_MO);
        else
            Beat(cntBeat).flag_n = 1;
        end
        
        if Beat(cntBeat).flag_n == false
            Beat(cntBeat).D = Beat(cntBeat).D1 + Beat(cntBeat).D2; % D = D1 + D2
        end
        
%         if Beat(cntBeat).flag_n == 0
%             locs_FP = locs_FP(locs_FP < Beat(cntBeat).idx_IRP - Beat(cntBeat).idx_S2 + rule.idx_thS2si);
%             if isempty(locs_FP) == false
%                 locs_FP = locs_FP(end) + Beat(cntBeat).idx_S2 - rule.idx_thS2si; % phase compantation (cuz extract data set)
%                 if (locs_FP > Beat(cntBeat).idx_AC) && (locs_FP < Beat(cntBeat).idx_IRP)
%                     Beat(cntBeat).D1 = Beat(cntBeat).seismo(Beat(cntBeat).idx_IRP) - ...
%                         Beat(cntBeat).seismo(locs_FP);
%                     Beat(cntBeat).D = Beat(cntBeat).D1 + Beat(cntBeat).D2;
%                 else
%                     Beat(cntBeat).flag_n = 1;
%                 end
%             end
%         end
    end
end
disp(['4) IRP,MO,AC detection : ',num2str(sum([Beat.flag_n])),'/',num2str(nBeat)]);

%% Part 6. Cal CITs
idx_En = ~[Beat.flag_n].*(1:length(Beat));
idx_En(idx_En == 0) = [];
for cntBeat = idx_En
    Beat(cntBeat).RRI = length(Beat(cntBeat).seismo) / fs;
    Beat(cntBeat).PEP = (Beat(cntBeat).idx_AO - Beat(cntBeat).idx_Q) / fs;
    Beat(cntBeat).ICT = (Beat(cntBeat).idx_AO - Beat(cntBeat).idx_MC) / fs;
    Beat(cntBeat).LVET = (Beat(cntBeat).idx_AC - Beat(cntBeat).idx_AO) / fs;
    Beat(cntBeat).IRT = (Beat(cntBeat).idx_MO - Beat(cntBeat).idx_AC) / fs;
end

%% Display Result
% Draw success or failure detection
figure;
for i = 1:nBeat
    v2 = [Beat(i).idx_raw(1)-Beat(1).idx_raw(1) -1; ...
        Beat(i).idx_raw(1)-Beat(1).idx_raw(1) 1; ...
        Beat(i).idx_raw(end)-Beat(1).idx_raw(1) 1; ...
        Beat(i).idx_raw(end)-Beat(1).idx_raw(1) -1];
    f2 = [1 2 3 4];
    if Beat(i).flag_n == 0
        patch('Faces',f2,'Vertices',v2,'FaceColor','green','FaceAlpha',.1);
    else
        patch('Faces',f2,'Vertices',v2,'FaceColor','red','FaceAlpha',.1);
    end
end

% Draw graph
hold on;
plot([Beat.seismo]./max([Beat.seismo])); 
plot([Beat.ecg]./max([Beat.ecg])*0.9);
plot([Beat.seismo_envelope]./max([Beat.seismo_envelope]),'--k');

% Draw FPs & text marker
for i = 1:nBeat
    plot(Beat(i).idx_MC +Beat(i).idx_raw(1)-Beat(1).idx_raw(1),Beat(i).seismo( Beat(i).idx_MC )./max([Beat.seismo]),'rv','MarkerFaceColor','r');
    plot(Beat(i).idx_ICP +Beat(i).idx_raw(1)-Beat(1).idx_raw(1),Beat(i).seismo( Beat(i).idx_ICP )./max([Beat.seismo]),'r^','MarkerFaceColor','r');
    plot(Beat(i).idx_AO +Beat(i).idx_raw(1)-Beat(1).idx_raw(1),Beat(i).seismo( Beat(i).idx_AO )./max([Beat.seismo]),'rv','MarkerFaceColor','r');
    plot(Beat(i).idx_AC +Beat(i).idx_raw(1)-Beat(1).idx_raw(1),Beat(i).seismo( Beat(i).idx_AC )./max([Beat.seismo]),'rv','MarkerFaceColor','r');
    plot(Beat(i).idx_IRP +Beat(i).idx_raw(1)-Beat(1).idx_raw(1),Beat(i).seismo( Beat(i).idx_IRP )./max([Beat.seismo]),'rv','MarkerFaceColor','r');
    plot(Beat(i).idx_MO +Beat(i).idx_raw(1)-Beat(1).idx_raw(1),Beat(i).seismo( Beat(i).idx_MO )./max([Beat.seismo]),'r^','MarkerFaceColor','r');
    plot(Beat(i).idx_R +Beat(i).idx_raw(1)-Beat(1).idx_raw(1),Beat(i).ecg( Beat(i).idx_R )./max([Beat.ecg])*0.9,'bv','MarkerFaceColor','b');
    plot(Beat(i).idx_Q +Beat(i).idx_raw(1)-Beat(1).idx_raw(1),Beat(i).ecg( Beat(i).idx_Q )./max([Beat.ecg])*0.9,'b^','MarkerFaceColor','b');
    
    text(Beat(i).idx_raw(1) - Beat(1).idx_raw(1), 1.04, ['#',num2str(i)],'fontsize',13);
    
    text(Beat(i).idx_MC +Beat(i).idx_raw(1)-Beat(1).idx_raw(1)-fs*0.03, Beat(i).seismo( Beat(i).idx_MC )./max([Beat.seismo])+0.045,'MC','fontsize',12);
    text(Beat(i).idx_ICP +Beat(i).idx_raw(1)-Beat(1).idx_raw(1)-fs*0.03, Beat(i).seismo( Beat(i).idx_ICP )./max([Beat.seismo])-0.04,'ICP','fontsize',12);
    text(Beat(i).idx_AO +Beat(i).idx_raw(1)-Beat(1).idx_raw(1)-fs*0.03, Beat(i).seismo( Beat(i).idx_AO )./max([Beat.seismo])+0.045,'AO','fontsize',12);
    text(Beat(i).idx_AC +Beat(i).idx_raw(1)-Beat(1).idx_raw(1)-fs*0.03, Beat(i).seismo( Beat(i).idx_AC )./max([Beat.seismo])+0.045,'AC','fontsize',12);
    text(Beat(i).idx_IRP +Beat(i).idx_raw(1)-Beat(1).idx_raw(1)-fs*0.03, Beat(i).seismo( Beat(i).idx_IRP )./max([Beat.seismo])+0.045,'IRP','fontsize',12);
    text(Beat(i).idx_MO +Beat(i).idx_raw(1)-Beat(1).idx_raw(1)-fs*0.03, Beat(i).seismo( Beat(i).idx_MO )./max([Beat.seismo])-0.04,'MO','fontsize',12);
    text(Beat(i).idx_R +Beat(i).idx_raw(1)-Beat(1).idx_raw(1)-fs*0.03, Beat(i).ecg( Beat(i).idx_R )./max([Beat.ecg])*0.9+0.045,'R','fontsize',12);
    text(Beat(i).idx_Q +Beat(i).idx_raw(1)-Beat(1).idx_raw(1)-fs*0.03, Beat(i).ecg( Beat(i).idx_Q )./max([Beat.ecg])*0.9-0.035,'Q','fontsize',12);
    text(Beat(i).idx_S1 +Beat(i).idx_raw(1)-Beat(1).idx_raw(1)+fs*0.03,Beat(i).seismo_envelope( Beat(i).idx_S1 )./max([Beat.seismo_envelope]),'S1','fontsize',12);
    text(Beat(i).idx_S2 +Beat(i).idx_raw(1)-Beat(1).idx_raw(1)+fs*0.03,Beat(i).seismo_envelope( Beat(i).idx_S2 )./max([Beat.seismo_envelope]),'S2','fontsize',12);
end
% xlim([77000 79000]);
% ylim([-1.5 1.5]);

% Draw CITs Result
figure;
subplot(5,1,1);
plot(locs_Rwave(idx_En)/fs/60,[Beat.RRI]*1000,'linewidth',1.5); 
y1 = ylabel('RRI (ms)'); set(y1,'rotation',45,'Position',[-0.2 mean(ylim) 0]);
set(gca,'XtickLabel',[]); set(gca,'xtick',[1:4]); xlim([-inf inf]); 

subplot(5,1,2);
plot(locs_Rwave(idx_En)/fs/60,[Beat.PEP]*1000,'linewidth',1.5); 
y2 = ylabel('PEP (ms)'); set(y2,'rotation',45,'Position',[-0.2 mean(ylim) 0]);
set(gca,'XtickLabel',[]); set(gca,'xtick',[1:4]); xlim([-inf inf]); 

subplot(5,1,3);
plot(locs_Rwave(idx_En)/fs/60,[Beat.ICT]*1000,'linewidth',1.5); 
y3 = ylabel('ICT (ms)'); set(y3,'rotation',45,'Position',[-0.2 mean(ylim) 0]);
set(gca,'XtickLabel',[]); set(gca,'xtick',[1:4]); xlim([-inf inf]); 

subplot(5,1,4);
plot(locs_Rwave(idx_En)/fs/60,[Beat.LVET]*1000,'linewidth',1.5); 
y4 = ylabel('LVET (ms)'); set(y4,'rotation',45,'Position',[-0.2 mean(ylim) 0]);
set(gca,'XtickLabel',[]); set(gca,'xtick',[1:4]); xlim([-inf inf]); 

subplot(5,1,5);
plot(locs_Rwave(idx_En)/fs/60,[Beat.IRT]*1000,'linewidth',1.5); 
y5 = ylabel('IRP (ms)'); set(y5,'rotation',45,'Position',[-0.2 mean(ylim) 0]);
set(gca,'xtick',[1:4]); xlim([-inf inf]); 
xlabel('Time (min)');

end