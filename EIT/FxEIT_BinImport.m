function [EIT_data] = FxEIT_BinImport(EIT_data_path,~,num_data)
if nargin > 1
    op = 'raw';
else
    op = 'mag';
end

scan_num = 512; % 1 raw data have 512 scan data
data_num = 256; % depend on protocol

cd(EIT_data_path);
dirlist = dir('.');

cnt_scan = 1;
for cnt_file = 1:length(dirlist)
    if length(dirlist(cnt_file).name) > 4
        if strcmp(dirlist(cnt_file).name(end-2:end),'bin')
            data_name{cnt_scan,1}=dirlist(cnt_file).name;
            cnt_scan = cnt_scan + 1;
        end
    end
end

% % time stamp
% for i = 1:length(data_name)
%     temp = data_name(i);
%     temp_time(i) = ConvText2Time(temp);
% end
% 
% subplot(211); plot(temp_time,'o');
% subplot(212); plot(diff(temp_time));
% drawnow;
% 
% temp = data_name(1);
% start_time = ConvText2Time(temp);
% temp = data_name(end);
% end_time = ConvText2Time(temp);
% if start_time > end_time
%     end_time = end_time + 24*60*60;
% end
% 
% time_index = linspace(start_time,end_time,(length(data_name)-1)*scan_num);
% inv_fs = mean(diff(time_index));
% fs = 1./inv_fs;
% temp = time_index(1)-inv_fs*scan_num:inv_fs:time_index(1)-inv_fs; % estimate first scan time info
% time_index = [temp time_index];
% clear dirlist;

if nargin < 3
    num_data = length(data_name);
end

EIT_data = zeros(data_num, num_data * scan_num);

cnt = 1;
for cnt_scan = 1:num_data-2
    file_path = strcat(EIT_data_path,'\',strcat('EIT_',num2str(cnt_scan),'.bin'));
%     file_path = strcat(EIT_data_path,'\',strcat()data_name{cnt_scan});
    fid = fopen(file_path);
    raw_data = fread(fid,data_num*scan_num*6,'uint8');
    raw_data = reshape(raw_data,data_num*6,scan_num);
    fclose(fid);
    
    for cnt_scan2 = 1:scan_num
        temp = raw_data(:,cnt_scan2);
        temp = reshape(temp,6,256);
        temp(1,:)=temp(1,:)-128;                     % 1. ch info
        temp(1:2,:) = temp([2 1],:);                 % 2. sat flag
        temp(3,:) = temp(4,:).*256 + temp(3,:);      % 3. R data
        temp(4,:) = temp(6,:).*256 + temp(5,:);      % 4. Q data
        temp(5:6,:) = []; % clear blank raw
        temp(temp>32767) = temp(temp>32767) - 65536; % conv 2's complement
        temp = temp';
        temp = temp([1:16:256 2:16:256 3:16:256 4:16:256 5:16:256 6:16:256 7:16:256 8:16:256 9:16:256 10:16:256 11:16:256 12:16:256 13:16:256 14:16:256 15:16:256 16:16:256],:);
        scan_data(:,:) = temp;

        if op == 'mag'
            EIT_data(:,cnt) = sqrt(scan_data(:,3).^2 + scan_data(:,4).^2);
        elseif op == 'raw'
            EIT_data(:,cnt) = scan_data(:,3) + 1j*scan_data(:,4);
        end
        clear temp;
        cnt = cnt + 1;
    end
    clear raw_data;
    
    disp(['data import : ' num2str(cnt_scan) ' / ' num2str(length(data_name))]);
end
    
% if num_data ~= length(EIT_data)
%     time_index(length(EIT_data)+1:end) = [];
% end
end

function [temp_time] = ConvText2Time(temp_text)
    tp = (strfind(temp_text,'_')); tp = tp{1};
    temp_time = 60*60*str2double(temp_text{1,1}(tp+1:tp+2)) + 60*str2double(temp_text{1,1}(tp+3:tp+4)) + str2double(temp_text{1,1}(tp+5:tp+6));

    % if cellfun('length',temp_text) < 26
    %     temp_time = 60*60*str2num(temp_text{1,1}(end-9:end-8)) + 60*str2num(temp_text{1,1}(end-7:end-6)) + str2num(temp_text{1,1}(end-5:end-4));
    % else
    %     temp_time = 60*60*str2num(temp_text{1,1}(end-12:end-11)) + 60*str2num(temp_text{1,1}(end-10:end-9)) + str2num(temp_text{1,1}(end-8:end-7));
    % end
end