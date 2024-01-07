function [idx,lack] = findidx(A,tp)
        if isdatetime(A(1))
            A = datenum(A);
        end
        if isdatetime(tp(1))
            tp = datenum(tp);
        end
        idx = NaN(size(tp));
        lack = NaN(size(tp));
        cnt_2 = 1;
        for cnt_1 = 2:length(A)
            if A(cnt_1) > tp(cnt_2)
                if A(cnt_1)-tp(cnt_2) > A(cnt_1-1)-tp(cnt_2)
                    idx(cnt_2) = cnt_1-1;
                    lack(cnt_2) = A(cnt_1-1)-tp(cnt_2);
                else
                    idx(cnt_2) = cnt_1;
                    lack(cnt_2) = A(cnt_1)-tp(cnt_2);
                end
                cnt_2 = cnt_2 + 1;
            end
            if cnt_2 > length(tp)
                idx(isnan(idx)) = length(A);
                break;
            end
        end
        idx(isnan(idx)) = length(A);
    end