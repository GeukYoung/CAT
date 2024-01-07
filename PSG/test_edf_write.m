function test_edf_write(filename,input_header,input_data)
n_channel = length(input_data); % channel number

% n_channel = 2;
% 
% %% input initial
% input_header.version = char(32*ones(1, 8));
% input_header.patientID = char(32*ones(1, 80));
% input_header.recordingID = char(32*ones(1, 80));
% input_header.startdate = char(32*ones(1, 8));
% input_header.starttime = char(32*ones(1, 8));
% input_header.nheaderrecord = zeros(1,0);
% input_header.records = zeros(1,0);
% input_header.duration = zeros(n_channel,1);
% input_header.labels = cell(n_channel,1);
% input_header.transducer = cell(n_channel,1);
% input_header.units = cell(n_channel,1);
% input_header.physMax = zeros(n_channel,1);
% input_header.physMin = zeros(n_channel,1);
% input_header.digMax = zeros(n_channel,1);
% input_header.digMin = zeros(n_channel,1);
% input_header.prefilt = cell(n_channel,1);
% input_header.n_sample = zeros(n_channel,1);
% input_header.samplerate = zeros(n_channel,1);
% 
% %%
% input_data = {0,0};
% input_data{1} = input_data1;
% input_data{2} = input_data2;
% 
% %% header information (input)
% filename = [filename '.edf']; % add '.edf'
% 
% % main info
% input_header.version = '0';
% input_header.patientID = 'X F X Female_33yr'; % [header.subject.ID ' ' header.subject.sex ' ' header.subject.birthdate ' ' header.subject.name]
% input_header.recordingID = 'Startdate 16-NOV-2003 X X X'; % ['Startdate ' num2str(header.day,'%02i') '-' monthstr '-' num2str(header.year,'%04i') ' ' header.ID ' ' header.technician ' ' header.equipment]);
% input_header.startdate = '16.11.03'; % date as dd.mm.yy
% input_header.starttime = '07.30.54'; % time as hh.mm.ss
% input_header.nheaderrecord = 256*(1+n_channel); % 2 : ch num
% input_header.duration = ones(n_channel,1)*1; % record duration (s)
% 
% % ch1 info
% input_header.labels{1} = 'EKG';
% input_header.transducer{1} = 'Ag/AgCl electrode';
% input_header.units{1} = 'mV';
% input_header.prefilt{1} = 'HP:0.5Hz LP:100Hz [enhanced cassette BW]';
% input_header.samplerate(1) = 200;
% 
% % ch2 info
% input_header.labels{2} = 'EIT';
% input_header.transducer{2} = 'Blue sensor';
% input_header.units{2} = 'AU';% 
% input_header.prefilt{2} = 'HP:0.1Hz LP:1Hz';
% input_header.samplerate(2) = 50;

%% initialize output header
output_header.version = char(32*ones(1, 8));
output_header.patientID = char(32*ones(1, 80));
output_header.recordingID = char(32*ones(1, 80));
output_header.startdate = char(32*ones(1, 8));
output_header.starttime = char(32*ones(1, 8));
output_header.nheaderrecord = char(32*ones(1, 8));
output_header.records = char(32*ones(1, 8));
output_header.duration = char(32*ones(1, 8));
output_header.channels = char(32*ones(1, 4));
output_header.labels = char(32*ones(n_channel, 16));
output_header.transducer = char(32*ones(n_channel, 80));
output_header.units = char(32*ones(n_channel, 8));
output_header.physMax = char(32*ones(n_channel, 8));
output_header.physMin = char(32*ones(n_channel, 8));
output_header.digMax = char(32*ones(n_channel, 8));
output_header.digMin = char(32*ones(n_channel, 8));
output_header.prefilt = char(32*ones(n_channel, 80));
output_header.n_sample = char(32*ones(n_channel, 8));

%% Data information (Record separation & scaling)
N = zeros(n_channel,1);
n_record = zeros(n_channel,1);
for i = 1:n_channel
    N(i) = length(input_data{i});
    n_record(i) = N(i)/input_header.samplerate(i);
    input_header.n_sample(i) = input_header.duration * input_header.samplerate(i);
end
records = max(n_record);

% Scale and convert to in16 (data)
maxV = zeros(n_channel,1); % phys
minV = zeros(n_channel,1); % phys

physMax = zeros(n_channel,1); % digi
physMin = zeros(n_channel,1); % digi
digMax = zeros(n_channel,1); % digi
digMin = zeros(n_channel,1); % digi
for i = 1:n_channel
    maxV(i) = max(input_data{i});
    minV(i) = min(input_data{i});
    % overflow check
    if maxV(i) > 3277 || minV(i) < -3277
        disp(strcat(num2str(i),' ch data > 16-bit integer, extreme values floored'));
        input_data{i} = input_data{i}./max([maxV(i) minV(i)])*3277;
%         input_data{i}(input_data{i} > 3277) = 3277; 
%         input_data{i}(input_data{i} < -3277) = -3277;
        maxV(i) = max(input_data{i});
        minV(i) = min(input_data{i});
    end
    maxdata = max([maxV(i) abs(minV(i))]);
    
    Scale = 32767 / maxdata;
    if maxdata < 1
        maxdata = 1;
    end
    
    physMax(i,:) = maxdata; 
    physMin(i,:) = -maxdata;
    digMax(i,:) = 32767;
    digMin(i,:) = -32767;
    input_data{i} = (input_data{i} .* Scale);
    
%     if ~strcmp(class(input_data{i}),'int16')
%         input_data{i} = int16(input_data{i});
%     end
end

DATAout = zeros(sum(input_header.samplerate.*input_header.duration),max(n_record));
Rs = ones(length(input_header.samplerate)+1,1);
Rs(2:end) = input_header.samplerate.*input_header.duration;
Rs=cumsum(Rs);
for i=1:length(input_data)
    temp_data = zeros(input_header.samplerate(i)*input_header.duration*max(n_record),1);
    temp_data(1:length(input_data{i})) = input_data{i};
    DATAout(Rs(i):(Rs(i+1)-1),:) = reshape(temp_data, input_header.samplerate(i).*input_header.duration, []);
end
DATAout = cast(DATAout,'int16');

if size(input_header.labels) ~= n_channel
  error 'Data dimension does not match header information';
end 

% Control for errors
if n_channel > 9999
  error 'Cannot write more than 9999 channels to an EDF file.';
end
if sum(N) > 99999999
  error 'Cannot write more than 99999999 data records (=samples) to an EDF file.';
end

%% make up header
%%% main header
% Version (8 ascii).
temp = '0';
output_header.version(1:length(temp)) = temp;

% Local patient (have been assured to be 80 ascii).
output_header.patientID(1:length(input_header.patientID)) = input_header.patientID;

% Recording identification (have been assured to be 80 ascii).
output_header.recordingID(1:length(input_header.recordingID)) = input_header.recordingID;

% Start date of recording (8 ascii).
output_header.startdate(1:length(input_header.startdate)) = input_header.startdate;

% Start time of recording (8 ascii)
output_header.starttime(1:length(input_header.starttime)) = input_header.starttime ;

% Number of bytes in header record (8 ascii); 256 for header + 256 for each data
temp = num2strc(input_header.nheaderrecord);
output_header.nheaderrecord(1:length(temp)) = temp;

% Number of data records (-1 if unknown; 8 ascii).
temp = num2strc(records);
output_header.records(1:length(temp)) = temp;

% Duration of a data record, in seconds (8 ascii).
temp = num2strc(input_header.duration(1));
output_header.duration(1:length(temp)) = temp;

% Number of signals (ns) in data record (4 ascii).
temp = num2strc(n_channel);
output_header.channels(1:length(temp)) = temp;

%%% channel header
% ns * Labels (have been assured to be 16 chars long).
for n = 1:n_channel
    if isempty(input_header.labels{n}) ~= 1
        ln = length(input_header.labels{n});
        if ln > 16    % label size over
            % label overflow
            ln = 16;
        end
        output_header.labels(n,1:ln) = input_header.labels{n}(1:ln);
    end
end

% ns * Transducer type (has been assured to be 80 chars).
for n = 1:n_channel
    if isempty(input_header.transducer{n}) ~= 1
        ln = length(input_header.transducer{n});
        if ln > 80    % label size over
            % label overflow
            ln = 80;
        end
        output_header.transducer(n,1:ln) = input_header.transducer{n}(1:ln);
    end
end

% ns * Create physdim-info (8 ascii).
for n = 1:n_channel
    if isempty(input_header.units{n}) ~= 1
        ln = length(input_header.units{n});
        if ln > 8    % label size over
            % label overflow
            ln = 8;
        end
        output_header.units(n,1:ln) = input_header.units{n}(1:ln);
    end
end

% % ns * Physical min/max & Digital min/max (8 ascii).
for i = 1:n_channel
    temp = num2strc(digMin(i));
    output_header.digMin(i,1:length(temp)) = temp;
    temp = num2strc(digMax(i));
    output_header.digMax(i,1:length(temp)) = temp;
    temp = num2strc(physMin(i));
    output_header.physMin(i,1:length(temp)) = temp;
    temp = num2strc(physMax(i));
    output_header.physMax(i,1:length(temp)) = temp;
end
% output_header.digMin = digMin;
% output_header.digMax = digMax;
% output_header.physMin = physMin;
% output_header.physMax = physMax;

%  ns * Prefiltering (80 ascii).
for n = 1:n_channel
    if isempty(input_header.prefilt{n}) ~= 1
        ln = length(input_header.prefilt{n});
        if ln > 80    % label size over
            % label overflow
            ln = 80;
        end
        output_header.prefilt(n,1:ln) = input_header.prefilt{n}(1:ln);
    end
end

%  ns * nr of samples in each data record (8 ascii).
for i = 1:n_channel
    temp = num2strc(input_header.n_sample(i));
    output_header.n_sample(i,1:length(temp)) = temp;
end

%% Write edf
% fid = fopen(filename); 
fid = fopen(filename, 'wb', 'ieee-le');

%%% >SUB: Write header record
% Version (1 byte).
fwrite(fid, output_header.version, 'char*1');

% Local patient (have been assured to be 80 chars).
fwrite(fid, output_header.patientID, 'char*1');

% Recording identification (have been assured to be 80 chars).
fwrite(fid, output_header.recordingID, 'char*1');

% Start date of recording (8 chars).
fwrite(fid, output_header.startdate, 'char*1');

% Start time of recording (8 chars).
fwrite(fid, output_header.starttime, 'char*1');

% Number of bytes in header record (8 chars); 256 for header + 256 for each data
fwrite(fid, output_header.nheaderrecord, 'char*1'); % number of bytes in header

% 1st reserved field (44 chars).
if isfield(input_header,'events')
    fprintf(fid, '%-44s', 'EDF+C'); % reserved (44 spaces)
else
    fprintf(fid, '%44s', ' '); % reserved (44 spaces)
end

% Number of data records (-1 if unknown; 8 chars).
fwrite(fid, output_header.records, 'char*1'); % number of records

% Duration of a data record, in seconds (8 chars).
fwrite(fid, output_header.duration, 'char*1'); % duration of record (=Fs)

% 4 ascii : number of signals (ns) in data record
fwrite(fid, output_header.channels, 'char*1'); % number of signals = channels

%%% >SUB: Write data record header
% Labels (have been assured to be 16 chars long).
fwrite(fid, output_header.labels', 'char*1'); % labels

% Transducer type (has been assured to be 80 chars).
fwrite(fid, output_header.transducer', 'char*1'); % transducer type (all spaces)

% Physical dimension (has been assured to be 8 chars).
fwrite(fid, output_header.units', 'char*1'); % phys dimension (all spaces)

% Physical min/max (assure 8 chars).
fwrite(fid, output_header.physMin', 'char*1'); % physical minimum
fwrite(fid, output_header.physMax', 'char*1'); % physical maximum

% Digital min/max (assure 8 chars).
fwrite(fid, output_header.digMin', 'char*1'); % digital minimum
fwrite(fid, output_header.digMax', 'char*1'); % digital maximum

% % Physical min/max (assure 8 chars).
% fprintf(fid, '%-8i', output_header.physMin);  % physical minimum
% fprintf(fid, '%-8i', output_header.physMax);  % physical maximum
% 
% % Digital min/max (assure 8 chars).
% fprintf(fid, '%-8i', output_header.digMin);  % digital minimum
% fprintf(fid, '%-8i', output_header.digMax);  % digital maximum

% ns * 80 ascii : ns * prefiltering (e.g. HP:0.1Hz LP:75Hz)
fwrite(fid, output_header.prefilt', 'char*1'); % prefiltering (all spaces)

% ns * 8 ascii : ns * nr of samples in each data record
fwrite(fid, output_header.n_sample', 'char*1'); % samples per record (= samplingrate)

% ns * 32 ascii : ns * reserved
fwrite(fid, 32*ones(32,n_channel), 'uint8'); % reserverd (32 spaces / channel)

%%% >SUB: Write data
fwrite(fid, DATAout, 'int16');
fclose(fid);
end

function str = num2strc(num) %#codegen
str = [];
flag_minus = 0;
if num < 0
    flag_minus = 1;
    num = -num;
end

while num > 0
    data = mod(num,10);
    str = [char(48+data), str]; %#ok
    num = (num-data)/10;
end

if flag_minus == 1
    str = [char(45), str];
end
end