function [Protocol, injPattern, meas] = func_generateEITDataMeasureProtocol(strDevice)

%% Voltage Measure Pattern (V+,V-)

meas = zeros(16,2);
meas(:,1) = 1:16;
meas(2:16,2) = 1:15;
meas(1,2) = 16;

%% Current Injection Pattern (Source, Sink)
if isequal(strDevice,'HV')
    chECG = [3 14]; % for HV
elseif isequal(strDevice,'AT')
    chECG = []; % for AT'
else
    Protocol = [];
    disp([strDevice ' is not supported'])
    return 
end

clear Protocol

cnt = 1;
for k=1:16
    
    if k == 1
        Sink = 16;
    else
        Sink = k-1;
    end
    
    if ~isempty(intersect(chECG,Sink))
        continue;
    end
    
    Source = Sink+1;
    if ~isempty(intersect(chECG,Source))
        Source = Source+1;
    end
    Source = mod(Source-1,16)+1;
    
    
    Protocol(cnt).injection = [Source Sink];
    Protocol(cnt).measure = meas;
    cnt = cnt+1;
end

if ~isempty(chECG)
    cnt = 15;
    Protocol(cnt).injection = [8 1];
    Protocol(cnt).measure = meas;
    
    cnt = 16;
    Protocol(cnt).injection = [16 9];
    Protocol(cnt).measure = meas;
end

tmp=cell2mat(struct2cell(Protocol));
injPattern = squeeze(tmp(1,:,:))';


