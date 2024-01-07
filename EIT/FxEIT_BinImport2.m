function [EIT_data] = FxEIT_BinImport2(EIT_data_path)
if nargin > 1
    op = 'raw';
else
    op = 'mag';
end

scan_num = 512; % 1 raw data have 512 scan data
data_num = 256; % depend on protocol

cd(EIT_data_path);
dirlist = dir('.');

flag.EIT = 1;
flag.HR = 1;
flag.SpO2 = 1;
cnt_EIT = 0;
cnt_HR = 0;
cnt_SpO2 = 0;
for cnt_file = 1:length(dirlist)
    if length(dirlist(cnt_file).name) > 4
        if strcmp(dirlist(cnt_file).name(end-3:end),'.bin')
            if strcmp(dirlist(cnt_file).name(1:4),'EIT_')
                cnt_EIT = cnt_EIT + 1;
                data_name_EIT{cnt_EIT,1}=dirlist(cnt_file).name;
            elseif strcmp(dirlist(cnt_file).name(1:3),'HR_')
                cnt_HR = cnt_HR + 1;
                data_name_HR{cnt_HR,1}=dirlist(cnt_file).name;
            elseif strcmp(dirlist(cnt_file).name(1:5),'SpO2_')
                cnt_SpO2 = cnt_SpO2 + 1;
                data_name_SpO2{cnt_SpO2,1}=dirlist(cnt_file).name;
            end
        end
    end
end
total_scan = cnt_EIT-1;

%% check exception
if (cnt_EIT == 0) || (cnt_EIT ~= cnt_HR) || (cnt_EIT ~= cnt_SpO2)
    if cnt_EIT == 0
        msgbox('data import fail');
        return;
    elseif cnt_HR == 0 && cnt_SpO2 == 0
        msgbox('only EIT data import');
        flag.HR = 0;
        flag.SpO2 = 0;
    elseif cnt_HR == 0
        msgbox('No HR data');
        flag.HR = 1;
    elseif cnt_SpO2 == 0
        msgbox('No SpO2 data');
        flag.SpO2 = 0;
    else
        msgbox('EIT, HR, SpO2 file number missmatch');
        flag.HR = 0;
        flag.SpO2 = 0;
    end
end


EIT_data.EIT_raw = zeros(data_num, total_scan * scan_num);

cnt = 1;
for cnt_EIT = 1:total_scan
    if flag.EIT == 1
        try
            file_path = strcat(EIT_data_path,'\EIT_',num2str(cnt_EIT-1),'.bin');
            fid = fopen(file_path);
            
            for cnt_scan2 = 1:scan_num
                EIT_data.version(:,cnt) = fread(fid,4,'char');
                EIT_data.Amp(:,cnt) = fread(fid,1,'double');
                EIT_data.Gain(:,cnt) = fread(fid,1,'double');
                EIT_data.Time(:,cnt) = fread(fid,1,'int64');
                EIT_data.Reserved(:,cnt) = fread(fid,44,'char');
                
                raw_data = fread(fid,data_num*6,'uint8');
                %         raw_data = reshape(raw_data,data_num*6,1);
                %         raw_data = fread(fid,data_num*scan_num*6,'uint8');
                %         raw_data = reshape(raw_data,data_num*6,scan_num);
                
                %         temp = raw_data(:,cnt_scan2);
                temp = reshape(raw_data,6,256);
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
                    EIT_data.EIT_raw(:,cnt) = sqrt(scan_data(:,3).^2 + scan_data(:,4).^2);
                elseif op == 'raw'
                    EIT_data.EIT_raw(:,cnt) = scan_data(:,3) + 1j*scan_data(:,4);
                end
                clear temp;
                cnt = cnt + 1;
                
            end
            clear raw_data;
            fclose(fid);
        catch
            fclose(fid);
            break;
        end
    end
    
    
    if flag.HR == 1
        file_path = strcat(EIT_data_path,'\HR_',num2str(cnt_EIT-1),'.bin');
        fid = fopen(file_path);
        raw_data = fread(fid,scan_num*10,'uint8');
        EIT_data.HR(scan_num*(cnt_EIT-1)+1:(scan_num*cnt_EIT)) = raw_data(1:10:end);
        fclose(fid);
    end
        
    if flag.SpO2 == 1
        file_path = strcat(EIT_data_path,'\SpO2_',num2str(cnt_EIT-1),'.bin');
        fid = fopen(file_path);
        raw_data = fread(fid,scan_num*9,'uint8');
        EIT_data.SpO2(scan_num*(cnt_EIT-1)+1:(scan_num*cnt_EIT)) = raw_data(1:9:end);
        fclose(fid);
    end
    
    disp(['data import : ' num2str(cnt_EIT) ' / ' num2str(total_scan)]);
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