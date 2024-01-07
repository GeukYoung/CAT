function [PV] = FxEIT_PVImport(PV_path)
fid = fopen(PV_path);
data = textscan(fid,'%s');

data{1,1}(1:29) = []; % remove header
data = reshape(data{1,1}(:),10,length(data{1,1})/10)';

ms = str2double(data(:,3));
PV.t_hms = datetime(str2num(data{1,1}(7:8))+2000 , str2num(data{1,1}(4:5)), str2num(data{1,1}(1:2)), ...
            str2num(data{1,2}(1:2)) , str2num(data{1,2}(4:5)), str2num(data{1,2}(7:8)), ms)';
        
% for i = 1:length(data)
%     PV.timeinfo(i) = 60*60*str2num(data{i,2}(1:2)) + 60*str2num(data{i,2}(4:5)) + str2num(data{i,2}(7:8));
% end
% 
% idx_prev = 1;
% for i = 2:length(data)
%     if PV.timeinfo(i) > PV.timeinfo(i-1)
%         PV.timeinfo(idx_prev:i) = linspace(PV.timeinfo(idx_prev),PV.timeinfo(i),i-idx_prev+1);
%         idx_prev = i;
%     end
% end

PV.t_hms = PV.t_hms(100:end);
PV.idx_Breath = str2double(data(100:end,4)); % Breath
PV.status = str2double(data(100:end,5)); % Status
PV.Paw = str2double(data(100:end,6)); % Paw
PV.Vaw = str2double(data(100:end,8)); % Volume

fclose(fid);
end
