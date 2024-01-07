function Y = catstruct(varargin)

narginchk(1,Inf) ;
N = nargin ;
if ~isstruct(varargin{end})
    if isequal(varargin{end},'sorted')
        narginchk(2,Inf) ;
        sorted = 1 ;
        N = N-1 ;
    else
        error('catstruct:InvalidArgument','Last argument should be a structure, or the string "sorted".') ;
    end
else
    sorted = 0 ;
end

FN = cell(N,1) ;
VAL = cell(N,1) ;
% parse the inputs
for cnt_s = 1:N
    X = varargin{cnt_s} ;
    if ~isstruct(X)
        error('catstruct:InvalidArgument',['Argument #' num2str(cnt_s) ' is not a structure.']) ;
    end
    
    FN{cnt_s} = fieldnames(X);
    VAL{cnt_s} = struct2cell(X);
end

for cnt_f = 1:length(VAL{1})
    if ~isstruct(VAL{1}{cnt_f})
        temp = [];
        for cnt_s = 1:N
            if size(VAL{cnt_s}{cnt_f},1) == 1
                temp = [temp; VAL{cnt_s}{cnt_f}];
            else
                temp = [temp, VAL{cnt_s}{cnt_f}];
            end
        end
        Y.(FN{1}{cnt_f}) = temp;
    else
        for cnt_s = 1:N
            VAL_iner{:,cnt_s} = struct2cell(VAL{cnt_s}{cnt_f});
            FN_iner{:,cnt_s} = fieldnames(VAL{cnt_s}{cnt_f});
        end
        for cnt_iners = 1:length(VAL_iner{:,1})
            temp = [];
            for cnt_s = 1:N
                if size(VAL_iner{cnt_s}{cnt_iners},1) < size(VAL_iner{cnt_s}{cnt_iners},2)
                    temp = [temp, VAL_iner{cnt_s}{cnt_iners}];
                else
                    temp = [temp; VAL_iner{cnt_s}{cnt_iners}];
                end
            end
            Y.(FN{1}{cnt_f}).(FN_iner{:,cnt_s}{cnt_iners}) = temp;
        end
        clear VAL_iner FN_iner;
    end
end