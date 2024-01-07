function [temp4] = autoYlim(varargin)
% narginchk(1, inf);

op_margin = 0.1;

data_obj = get(gca,'Children');
data = get(data_obj,'YData');
if iscell(data)
    data = [data{:}];
end

numOrigInputArgs = numel(varargin);
if numOrigInputArgs ~= 0
    while numOrigInputArgs
        switch varargin{1}
            case 'min'
                op_ymin = varargin{2};
                varargin(1:2) = [];
                numOrigInputArgs = length(varargin);
            case 'max'
                op_ymax = varargin{2};
                varargin(1:2) = [];
                numOrigInputArgs = length(varargin);
            case 'range'
                if length(varargin{2}) == 2
                    op_range = varargin{2};
                    data = data(op_range(1):op_range(2));
                else
                    op_range = varargin{2};
                    data = data(op_range:end);
                end
                varargin(1:2) = [];
                numOrigInputArgs = length(varargin);
            case 'margin'
                if length(varargin{2}) == 2
                    temp_lim = varargin{2};
                    y_margin_min = temp_lim(1);
                    y_margin_max = temp_lim(2);
                else
                    op_margin = varargin{2};
                end
                varargin(1:2) = [];
                numOrigInputArgs = length(varargin);
            otherwise
                msgbox(['Unvalid input :' varargin{1}]);
                numOrigInputArgs = 0;
        end
    end
end
ylim_origin = [min(data) max(data)];

if ~exist('y_margin_min')
    y_margin_min = op_margin*diff(ylim_origin);
    y_margin_max = op_margin*diff(ylim_origin);
end
temp3 = [ylim_origin(1)-y_margin_min ylim_origin(2)+y_margin_max];

%% 
temp4 = [round(temp3(1),-floor(log10(abs(temp3(1))))+2) ...
    round(temp3(2),-floor(log10(abs(temp3(2))))+2)];

if exist('op_ymin')
    temp4(1) = op_ymin;
end
if exist('op_ymax')
    temp4(2) = op_ymax;
end
    
ylim(temp4);
yticks(temp4);

end
