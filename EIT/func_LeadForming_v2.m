function [vLeadForming, indxRef] = func_LeadForming_v2(filt, filtHigh, RpeakFlag)


%% Use case

% tp_lead = [1 40*100];
% 
% indRpeak_leadforming = qrs_i_raw( qrs_i_raw>= tp_lead(1) & qrs_i_raw<= tp_lead(2));
% 
% RpeakFlag = false(indRpeak_leadforming(end),1);
% RpeakFlag(indRpeak_leadforming) = true;
% 
%
% fs = 100;
% N = 2; fc = 8; butter = 1; bpf = 1;
% if 1
%     filt = FxEIT_Filter(ohm_lead,fs,N,fc,butter,bpf);
% end
% 
% N = 2; fc = 0.5; butter = 1; bpf = 2;
% if 1
%     filtHigh = FxEIT_Filter(filt,fs,N,fc,butter,bpf);
% end

% [vLeadForming, indxRef] = func_LeadForming_v2(filt,filtHigh, RpeakFlag);
% 
% CVS2 = vLeadForming'*filtHigh;


%%
indxRef = find(RpeakFlag == 1);


indCycle = zeros(length(indxRef)-1,2);
indCycle(:,1) = indxRef(1:end-1);
indCycle(:,2) = indxRef(2:end);



isDiffLeadForming = true;

indxRef_Validated =unique(indCycle(:));

r= 0.50;
indCycleES = round(indCycle(:,1)*r+indCycle(:,2)*(1-r));


if isDiffLeadForming
    data = filt;
    
    r= 0.50;
    dataFidel = data(:,indCycle(:,1))-data(:,indCycleES); % ED, ES
    data = filtHigh;
else
    
    data = filtHigh;
    dataFidel = data(:,indxRef_Validated);
end


data = filt;

dataDiff = [];

for i=1:size(indCycle,1)-2
    
    ind = indCycle(i,1):indCycle(i,2);
    ind2 = indCycle(i+1,1):indCycle(i+1,2);
    
    
    if length(ind) > length(ind2)
        indInterp = ind2;
        
    else
        indInterp = ind;
    end
    
    if length(indInterp)==1
        continue;
    end
    
    indInterp = indInterp - indInterp(1) + 1;
    
    dataTmp1 = zeros( size(data,1), length(indInterp));
    dataTmp2 = zeros( size(data,1), length(indInterp));
    
    for k=1:size(data,1)
        
        if length(ind) > length(ind2)
            dataTmp1(k,:) = interp1((ind-ind(1))/(ind(end)-ind(1)), data(k,ind), (indInterp-indInterp(1))/(indInterp(end)-indInterp(1)));
            dataTmp2(k,:) = data(k,ind2);
        else
            dataTmp1(k,:) = data(k,ind);
            dataTmp2(k,:) = interp1((ind2-ind2(1))/(ind2(end)-ind2(1)), data(k,ind2), (indInterp-indInterp(1))/(indInterp(end)-indInterp(1)));
        end
    end
    
    dataDiffTmp = (dataTmp1 - dataTmp2);
    dataDiffTmp(:,[1 end]) = [];
    
    
    dataDiff = [dataDiff dataDiffTmp];
    
end





EIT = dataFidel;
EIT2 = dataDiff;

isUseHalfData = false;

% [vLeadForming] = func_LeadForming(dataFidel,dataDiff, false);
%%
% isUseHalfData = true;

if nargin == 1
    EIT2 = [];
    isUseHalfData = true;
end

% EIT : 208 times N

Nchannel = 16;
%% find reciprocal pairs
inj = repmat((1:Nchannel)',1, Nchannel )';
mea = inj';

inj = inj(:);
mea = mea(:);

rmv_indx=func_rmv_skip(Nchannel,0);

inj(rmv_indx)=[];
mea(rmv_indx)=[];


indpair=sort([inj mea],2);

indReci = zeros(size(indpair,1),2);

for i = 1:length(indReci)
    indReci(i,:) = find( indpair(i,1) == indpair(:,1) & indpair(i,2) == indpair(:,2));
end

indChReciAll = unique(sort(indReci,2),'rows');


% r= 0.50;
% indxRef2=round((indxRef(1:end-1)*r+indxRef(2:end)*(1-r)));

if isUseHalfData
    dataFidel = EIT(indChReciAll(:,1),:) + EIT(indChReciAll(:,2),:);
else
    dataFidel = EIT;
end

% if 1
%     dataFidel = dataFidel(:,indxRef(1:end-1))-dataFidel(:,indxRef2); % ED, ES
% end


if ~isempty(EIT2)
    if isUseHalfData
        dataFidel2 = EIT2(indChReciAll(:,1),:) + EIT2(indChReciAll(:,2),:);
    else
        dataFidel2 = EIT2;
    end
    B = dataFidel2';
end

A=dataFidel';

if isempty(EIT2)
    x = A\ones(size(A,1),1);
else
    x = [A; B]\[ones(size(A,1),1) ; zeros(size(B,1),1)];
end


if isUseHalfData
    vLeadForming = zeros(size(EIT,1),1);
    vLeadForming(indChReciAll(:,1)) = x;
    vLeadForming(indChReciAll(:,2)) = x;
else
    vLeadForming = x;
end


vLeadForming = vLeadForming/norm(vLeadForming);
