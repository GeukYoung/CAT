function [Biopac] = FxImport_Biopac(path_Biopac)
    op.ch_Trig = 1;
    op.ch_ECG = 3;
    
    fid = fopen(path_Biopac);
    temp_txt = textscan(fid,'%f %f %f');
    fclose(fid);
    
    if (std(temp_txt{op.ch_Trig}(1:2:end)) > std(temp_txt{op.ch_Trig}(2:2:end))) || isnan(std(temp_txt{op.ch_Trig}(2:2:end)))
        Biopac.Trig = temp_txt{op.ch_Trig}(1:2:end);
        Biopac.ECG = temp_txt{op.ch_ECG}(1:2:end);
    end
    op.n_length = length(Biopac.Trig);
    temp = abs(fft(Biopac.Trig));
    temp = temp(1:op.n_length*0.5);
    temp(1:100) = 0; % rem DC
    figure; plot(temp);
    [~,tp] = max(temp);
    fs_1 = round(op.n_length / tp / 2); % find 1/fs
    
    temp = movmax(Biopac.Trig,[fs_1 fs_1]);
%     figure; plot(Biopac.Trig); hold on; plot(temp);
    
    th = std(temp);
    temp = diff(temp);
    
    tp_max = find(temp>th);
    tp_min = find(temp(tp_max(1):end)<-th) + tp_max(1);
    if length(tp_min) < length(tp_max)
        tp_min(end+1) = length(temp);
    end
    tp_length = tp_min - tp_max;
    [~,tp_sEIT] = max(tp_length);
    idx_sEIT = tp_max(tp_sEIT)-fs_1*10;
    
    figure;
    subplot(211); plot(Biopac.Trig); hold on; bar(idx_sEIT,max(Biopac.Trig),'r');
    subplot(212); plot(Biopac.Trig); hold on; bar(idx_sEIT,max(Biopac.Trig),'r'); xlim([idx_sEIT-fs_1*100 idx_sEIT+fs_1*100])
    
    Biopac.Trig(1:idx_sEIT) = [];
    Biopac.ECG(1:idx_sEIT) = [];
end

