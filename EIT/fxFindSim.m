function [simdata] = fxFindSim(data)
fmdl = ng_mk_cyl_models([0.5,1,0.05],[16,0.25],[0.05]); 
fmdl.electrode = fmdl.electrode([6:16 , 1:5]);
% show_fem(fmdl,[0 1]);
fmdl.stimulation = mk_stim_patterns(16,1,[0,1],[0,1],{},1);

cnt = 1;
for i = 1:10:100
    img = mk_image(fmdl,i);
    temp = fwd_solve(img);
    vh(:,cnt) = temp.meas/min(temp.meas);
    % plot(vh)
    % plot(mean(data,2)/max(mean(data,2)));
    xc(cnt) = max(xcorr(vh(:,cnt),mean(data,2))/min(mean(data,2)));
    disp(num2str(cnt))
    cnt = cnt + 1;
end
[~,index] = max(xc);
simdata = vh(:,index);
plot(mean(data,2)/min(mean(data,2)));
hold on; plot(simdata,'r');
plot(vh);