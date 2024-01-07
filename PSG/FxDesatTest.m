function [Desat_flag,DesatEvent] = FxDesatTest(PSG,opt)
% input
%   PSG : SPO2 data among PSG data set
%   opt : draw result option
% output
%   Desat_flag : normal->0, desat->1
%   DesatEvent

    th_Desat = 4; % 4% drop
    tolerance_Desat = 1; % 1% tolerance
    SPO2 = round(PSG.filt_data);
    
    event_start = 0;
    event_end = 0;
    state = 1;
    SPO2_Variation = 0;
    cntEvent = 1;
    for cnt = 2:length(SPO2)
        switch state
            case 1 % Resting state
                if SPO2(cnt) < SPO2(cnt-1) % SPO2 drop
                    event_start = cnt; % Start time storing
                    event_end = cnt;
                    state = 2;
                    SPO2_Variation = 0;
                end
            case 2 % Event state
                if SPO2(cnt) < SPO2(cnt-1) % SPO2 drop
                    SPO2_Variation = 0;
                elseif SPO2(cnt) > SPO2(cnt-1) % SPO2 increase
                    if SPO2_Variation >= tolerance_Desat
                        if SPO2(event_start-1) - SPO2(event_end-1) >= th_Desat
                           DesatEvent(cntEvent,:) = [event_start event_end];
                           cntEvent = cntEvent + 1;
                        end
                        SPO2_Variation = 0;
                        state = 1;
                    else
                        event_end = cnt - 1; % End time storing
                        SPO2_Variation = SPO2_Variation + 1;
                    end
                end
        end
    end

Desat_flag = zeros(size(SPO2));
for cnt = 1:length(DesatEvent)
    Desat_flag(DesatEvent(cnt,1):DesatEvent(cnt,2)) = 1;
end

% draw result

if nargin > 1
    figure('name','SpO2 Desaturation'); set(gcf,'units','normalized','outerposition',[0 .3 1 0.4]);
    hold on;
    for cnt = 1:length(DesatEvent)
        vert = [DesatEvent(cnt,1) 80; ...
            DesatEvent(cnt,1) 100; ...
            DesatEvent(cnt,2) 100; ...
            DesatEvent(cnt,2) 80];
        fac = [1 2 3 4];
        patch('Faces',fac,'Vertices',vert,'FaceColor','red','FaceAlpha',.7);
%         colormap jet; caxis([-1.5 1.5]);
    end
    plot(SPO2,'linewidth',1.5); ylim([80 100]); xlim([-inf inf]);
    set(gcf, 'color', 'white');
    ylabel('SpO2 (%)','fontsize',14);
end
end
    
