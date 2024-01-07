function imdl = imdl_human(type_elec)
if nargin < 1
    type_elec = 360;
else
    if type_elec ~= 360
        type_elec = 220;
    end
end

% clear all;
% EIDORS_path = 'D:\OneDrive\1. Data\99. EIT SW\EIDORS\eidors-v3.8_2dGR\eidors\startup.m';
% EITtool_path = 'D:\OneDrive\1. Data\99. EIT SW\EIT tool\EITtools.m';
% 
% run(EIDORS_path);
% run(EITtool_path);

elec = [105.211491442543,9.61002444987770;54.7762836185819,23.6051344743276;19.3924205378973,54.7640586797066;8.30195599022004,93.3166259168704;8.30195599022004,132.397310513447;24.9376528117359,171.477995110024;64.8105134474327,198.411980440098;116.830073349633,208.974327628362;168.849633251834,208.710268948655;220.869193154034,198.676039119804;260.477995110024,171.477995110024;277.641809290953,132.925427872861;277.641809290953,93.5806845965770;266.551344743276,55.0281173594132;231.167481662592,24.1332518337408;180.996332518337,9.87408312958431];
elec(:,2) = -elec(:,2);
elec = elec([9:16 1:8],:);
figure;
axis image; hold on;
for i = 1:16
    plot(elec(i,1),elec(i,2),'ro','MarkerSize',10,'MarkerFaceColor','r');
    text(elec(i,1),elec(i,2),num2str(i),'color','b','FontWeight','bold','FontSize',14)
end
op.hmax = 0.03;
op.thick = 0.01;
imdl = FxEIT_FER2(elec,op);

figure; show_fem(imdl.fwd_model,[0 1]);

if type_elec == 220
    %% find electrode position
    % find 12clk point
    [~, bnd_idx_new] = min(pdist2(imdl.fwd_model.nodes, ...
        [(max(imdl.fwd_model.nodes(:,1))+min(imdl.fwd_model.nodes(:,1)))/2 max(imdl.fwd_model.nodes(:,2))]));
    figure; show_fem(imdl.fwd_model,[0 1]); hold on; plot(imdl.fwd_model.nodes(bnd_idx_new,1),imdl.fwd_model.nodes(bnd_idx_new,2),'ro','markerfacecolor','r');
    
    % find 2nd bnd point (clock wise dir)
    bnd_idx_0 = imdl.fwd_model.boundary;
    temp = find(bnd_idx_0==bnd_idx_new,2);
    for i = 1:2
        if temp(i) < size(bnd_idx_0,1)
            bnd_idx_next(i,1) = bnd_idx_0(temp(i),2);
            bnd_idx_next(i,2) = temp(i);
        else
            bnd_idx_next(i,1) = bnd_idx_0((temp(i)-size(bnd_idx_0,1)),1);
            bnd_idx_next(i,2) = (temp(i)-size(bnd_idx_0,1));
        end
    end
    if imdl.fwd_model.nodes(bnd_idx_next(1,1),1) > imdl.fwd_model.nodes(bnd_idx_next(2,1),1) % clock wise
        bnd_idx_new(2) = bnd_idx_next(1,1);
        bnd_idx_0(bnd_idx_next(1,2),:) = 0;
    else
        bnd_idx_new(2) = bnd_idx_next(2,1);
    end
    hold on; plot(imdl.fwd_model.nodes(bnd_idx_new(2),1),imdl.fwd_model.nodes(bnd_idx_new(2),2),'ro','markerfacecolor','r');
    
    % find others
    for i = 3:length(bnd_idx_0)+1
        temp = find(bnd_idx_0==bnd_idx_new(i-1),1);
        if temp < size(bnd_idx_0,1)
            bnd_idx_new(i) = bnd_idx_0(temp,2);
            bnd_idx_0(temp,:) = 0;
        else
            bnd_idx_new(i) = bnd_idx_0((temp-size(bnd_idx_0,1)),1);
            bnd_idx_0((temp-size(bnd_idx_0,1)),:) = 0;
        end
        plot(imdl.fwd_model.nodes(bnd_idx_new(i),1),imdl.fwd_model.nodes(bnd_idx_new(i),2),'ro','markerfacecolor','r');
        pause(0.001);
    end
    
    bnd_pnt = imdl.fwd_model.nodes(bnd_idx_new,:);
    bnd_dist = dist([bnd_pnt; bnd_pnt(1,:)]');
    bnd_dist = diag(bnd_dist,1);
    bnd_dist_sum = sum(bnd_dist);
    bnd_dist_cumsum = cumsum(bnd_dist);
    
    %% fwd_model
    dim.c2e = 3; % 중심부터 첫전극
    dim.e2e = 4; % 전극간 간격
    dim.l2l = 4; % 전극 패드사이 여유길이 (< 11 이하)
    dim.ldomain = 105; % 둘레
    
    elec_pos = ([(0:dim.e2e:dim.e2e*3) ((0:dim.e2e:dim.e2e*3)+dim.l2l+dim.e2e*3)]+dim.c2e)/dim.ldomain*bnd_dist_sum;
    % elec_zdist = dim.l2l/dim.ldomain*bnd_dist_sum;
    % elec_r = 1.5/2/dim.ldomain*bnd_dist_sum;
    for i = 1:length(elec_pos)
        [~, elec_idx(i)] = min(abs(bnd_dist_cumsum-elec_pos(i)));
    end
    for i = 1:length(elec_pos)
        [~, elec_idx(i+8)] = min(abs(bnd_dist_cumsum-(bnd_dist_sum-elec_pos(i))));
    end
    
    % distance round up down error calibration....
    elec_idx = elec_idx+1;
    elec_idx = elec_idx([8:-1:1 9:16]);
    figure;
    show_fem(imdl.fwd_model); hold on; plot(bnd_pnt(elec_idx,1),bnd_pnt(elec_idx,2),'ro','MarkerFaceColor','r');
    for i = 1:16
        text(bnd_pnt(elec_idx(i),1),bnd_pnt(elec_idx(i),2),num2str(i),'color','b','FontWeight','bold','FontSize',14)
    end
    axis off image;
    
    op.bnd_pos = bnd_pnt(1:2:end,:);
    op.bd_n = 1;
    op.bd_alpha=1e-11;
    op.thick=0.025;
    
    imdl = FxEIT_FER2(bnd_pnt(elec_idx,:),op);
end
end

