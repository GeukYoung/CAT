function [EIT_Voltage] = FxATRM100_VoltageImport(EIT_data_path)
meas_num = 256; % 256 measurement

cd(EIT_data_path);
cnt = 1;

file_path = strcat(EIT_data_path,'\','VoltageGain_0.bin');
fid = fopen(file_path);
raw_data = fread(fid,'double');
scan_num = size(raw_data,1)/(meas_num);
raw_data = reshape(raw_data,meas_num,scan_num);

fclose(fid);
EIT_Voltage = raw_data;

clear raw_data;