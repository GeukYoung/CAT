function [Biopac] = FxImport_Biopac_old(path_Biopac)
    op.ch_Trig = 2;
%     op.ch_Transonic = 2;
    op.ch_ECG = 3;

    fid = fopen(path_Biopac);
    temp_txt = textscan(fid,'%f %f %f');
    fclose(fid);
    figure; plot(temp_txt{op.ch_Trig});
    if sum(temp_txt{op.ch_Trig}) < 0
        temp_txt{op.ch_Trig} = -temp_txt{op.ch_Trig};
    end
    
    if (std(temp_txt{op.ch_Trig}(1:2:end)) > std(temp_txt{op.ch_Trig}(2:2:end))) || isnan(std(temp_txt{op.ch_Trig}(2:2:end)))
        Biopac.Trig = temp_txt{op.ch_Trig}(1:2:end);
%         Biopac.Transonic = temp_txt{op.ch_Transonic}(1:2:end);
        Biopac.ECG = temp_txt{op.ch_ECG}(1:2:end);
    else
        Biopac.Trig = temp_txt{op.ch_Trig}(2:2:end);
%         Biopac.Transonic = temp_txt{op.ch_Transonic}(2:2:end);
        Biopac.ECG = temp_txt{op.ch_ECG}(2:2:end);
    end
    
    if mean(Biopac.Trig(1:100)) > 1.5
        Biopac.Trig(1:find(Biopac.Trig<1.5,1,'first')) = 0;
    end
%     figure;
%     plot(Biopac.Trig);
%     
%     figure;
%     plot(abs(diff(Biopac.Trig)));
    
    temp = abs(diff(Biopac.Trig));
    
    op.n_length = length(Biopac.Trig);
    temp = abs(fft(Biopac.Trig));
    temp = temp(1:op.n_length*0.5);
    temp(1:100) = 0; % rem DC
    figure; plot(temp);
    [~,tp] = max(temp);
    fs_1 = round(op.n_length / tp / 2); % find 1/fs
    
    temp = movmax(Biopac.Trig,[fs_1 fs_1]);
    temp(temp<1.5) = 0;
    temp(temp>1.5) = 1;
    
    temp_edge = diff(temp);
    temp_edge_u = temp_edge;
    temp_edge_u(temp_edge_u<0) = 0;
    [temp_edge_u] = find(temp_edge_u);
    temp_edge_d = -temp_edge;
    temp_edge_d(temp_edge_d<0) = 0;
    [temp_edge_d] = find(temp_edge_d);
    if size(temp_edge_d,1) > 0
        if temp_edge_d(1) < temp_edge_u(1)
            temp_edge_d(1) = [];
        end
    end
    if length(temp_edge_d) ~= length(temp_edge_u)
        temp_edge_d = [temp_edge_d; length(temp)];
    end
    figure; plot(temp); hold on;
    plot(temp_edge_u,temp(temp_edge_u),'r^'); plot(temp_edge_d,temp(temp_edge_d),'bv'); ylim([-0.1 1.1])

    temp_length = temp_edge_d - temp_edge_u;
    [~, tp] = max(temp_length);
%     figure; plot(Biopac.Trig); hold on; plot(temp);
    idx_sEIT = temp_edge_u(tp) - fs_1 * 2;
%     th = std(temp);
%     temp = diff(temp);
%     
%     tp_max = find(temp>th);
%     tp_min = find(temp(tp_max(1):end)<-th) + tp_max(1);
%     if length(tp_min) < length(tp_max)
%         tp_min(end+1) = length(temp);
%     end
%     tp_length = tp_min - tp_max;
%     [~,tp_sEIT] = max(tp_length);
%     idx_sEIT = tp_max(tp_sEIT)-fs_1*10;
    
    figure;
    subplot(211); plot(Biopac.Trig); hold on; bar(idx_sEIT,max(Biopac.Trig),'r'); plot(temp);
    subplot(212); plot(Biopac.Trig); hold on; bar(idx_sEIT,max(Biopac.Trig),'r'); xlim([idx_sEIT-fs_1*100 idx_sEIT+fs_1*100])
    
    Biopac.Trig(1:idx_sEIT) = [];
%     Biopac.Transonic(1:idx_sEIT) = [];
    Biopac.ECG(1:idx_sEIT) = [];
    Biopac.Fs = fs_1*100;
end

