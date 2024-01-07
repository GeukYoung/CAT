function [chMask, chRemain, chPairReci, chReci_unique, injmeaTable, Protocol, injPattern, measPattern] = func_get_channel_Mask_Reciprocity(strDevice, isChannelNumberingTo256)

if nargin == 1
    isChannelNumberingTo256 = [];
end


if isempty(isChannelNumberingTo256)
   isChannelNumberingTo256 = false; 
end


[Protocol, injPattern, measPattern] = func_generateEITDataMeasureProtocol(strDevice);

if isempty(Protocol)
    disp([strDevice ' is not supported'])
    return
end


tmp=cell2mat(struct2cell(Protocol));
injPattern = squeeze(tmp(1,:,:))';


%% Masking

injmeaTable = [];

for k=1:length(Protocol)
    mea = Protocol(k).measure;
    inj = Protocol(k).injection;
    injmeaTable=[injmeaTable; [repmat(inj,length(mea),1) mea]];
end

rmvflag = false(size(injmeaTable,1),1);

for k=1:length(rmvflag)
    if injmeaTable(k,1) == injmeaTable(k,3) || injmeaTable(k,1) == injmeaTable(k,4)
        rmvflag(k) = true;
    end
    
    if injmeaTable(k,2) == injmeaTable(k,3) || injmeaTable(k,2) == injmeaTable(k,4)
        rmvflag(k) = true;
    end
    
end

chMask = find(rmvflag);
chRemain = find(~rmvflag);

%% Reciprocity Channel Pairs

% isChannelNumberingTo256 = false;

if ~isChannelNumberingTo256
    injmeaTable(chMask,:)=[];
end
channel = (1:size(injmeaTable,1))';

chReci = zeros(size(injmeaTable,1),1);
for k=1:size(injmeaTable,1)
    
    ref1 = injmeaTable(k,1:2);
    ref2 = injmeaTable(k,3:4);
    [rr, ~]=find(ref1(1) == injmeaTable(:,3) & ref1(2) == injmeaTable(:,4) & ref2(1) == injmeaTable(:,1) & ref2(2) == injmeaTable(:,2));
    rr = unique(rr);
    if isempty(rr) || rr == k
        chReci(k) = nan;
    else
        chReci(k) = rr;
    end
    
    if isChannelNumberingTo256 && rmvflag(k)
        chReci(k) = nan;
    end
    
end


indReci = find(~isnan(chReci));

chPairReci= [indReci chReci(indReci)];

chReci_unique = unique(sort(chPairReci,2),'rows');
%% Save
if 0
    if isempty(chECG)
        save Protocol_AirTom.mat Protocol injPattern chMask chRemain chPairReci
    else
        save Protocol_HemoVista.mat Protocol injPattern chMask chRemain chPairReci
    end
end

