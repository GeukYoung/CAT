function [vLeadForming] = func_LeadForming(EIT, EIT2, isUseHalfData)


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