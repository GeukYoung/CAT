function [Step,Section] = FxParam_findPEEP(PEEP)

opt.th_section = 0.8; % threshold of PEEP section
opt.nTV_section = 3; % minimum TV number of PEEP section
opt.th_step = 1; % threshold of PEEP step
opt.min_nstep = 3; % minimum PEEP number of PEEP step

cnt.section = 1;
cnt.step = 1;
flag.update_section = 0;
flag.update_step = 1;

Section.idx = [];
Section.PEEP = [];
Section.dir = [];
Section.Pdiff = [];

Step(1).idx = [];
Step(1).dir = [];

for cnt_TV = 1:length(PEEP)
    if cnt_TV >= opt.min_nstep
        % calc PEEP difference between N & N-nTV_section+1
        diff_PEEP = abs(PEEP(cnt_TV)-PEEP(cnt_TV-opt.nTV_section+1));
        
        % find PEEP section
        if diff_PEEP < opt.th_section % check PEEP level change over th lvl
            % check section overlap
            if cnt.section > 1 % atleast 2nd PEEP section needed
                if cnt_TV - Section.idx(cnt.section-1) < 2
                    cnt.section = cnt.section - 1; % give same cnt.section by -1
                    flag.update_section = 0;
                else
                    flag.update_section = 1; % new section start info
                end
            end
            
            % stored PEEP lv & idx
            Section.idx(cnt.section) = cnt_TV;
            Section.PEEP(cnt.section) = PEEP(cnt_TV);
            
            % stored PEEP direction
            if cnt.section > 1 % atleast 2nd PEEP section needed
                if Section.PEEP(cnt.section) > Section.PEEP(cnt.section-1)
                    Section.dir(cnt.section) = 1; % increase
                    Section.Pdiff(cnt.section) = abs(Section.PEEP(cnt.section) - Section.PEEP(cnt.section-1));
                else
                    Section.dir(cnt.section) = -1; % decrease
                    Section.Pdiff(cnt.section) = abs(Section.PEEP(cnt.section) - Section.PEEP(cnt.section-1));
                end
            end
            
            % find PEEP step
            if flag.update_section == 1 % check section update
                if (Section.Pdiff(cnt.section) > opt.th_step) && (Section.Pdiff(cnt.section-1) > opt.th_step) && (Section.dir(cnt.section)==Section.dir(cnt.section-1)) % check same direction (increase or decrease)
                    if flag.update_step % first detection (1 & 2 & 3)
                        Step(cnt.step).idx = [Section.idx(cnt.section-2) Section.idx(cnt.section-1) Section.idx(cnt.section)];
                        Step(cnt.step).dir = Section.dir(cnt.section);
                        flag.update_step = 0;
                    else % (4 & ...)
                        Step(cnt.step).idx(end) = []; % update cnt.section-1
                        Step(cnt.step).idx = [Step(cnt.step).idx Section.idx(cnt.section-1) Section.idx(cnt.section)];
                    end
                else % end of step
                    if flag.update_step == 0
                        if length(Step(cnt.step).idx) < opt.min_nstep
                            Step(cnt.step).idx = []; % clear stored PEEP step
                            Step(cnt.step).dir = []; % clear stored PEEP step
                        else
                            Step(cnt.step).idx = Step(cnt.step).idx - 1;
                            cnt.step = cnt.step + 1;
                        end
                        flag.update_step = 1; % reset first detection flag
                    end
                end
            end
            
            cnt.section = cnt.section + 1;
        end
    end
end

%%

figure;
tiledlayout(5, 1);
h(1) = nexttile;
plot(PEEP,'.'); hold on;
plot(Section.idx,PEEP(Section.idx),'r.'); ylabel('PEEP data');
h(2) = nexttile;
plot(Section.idx,Section.idx,'.'); ylabel('section idx')
h(3) = nexttile;
plot(Section.idx,Section.PEEP,'.'); ylabel('section PEEP')
h(4) = nexttile;
plot(Section.idx,Section.dir,'.'); ylabel('section dir')
h(5) = nexttile;
plot(Section.idx,Section.Pdiff,'.'); ylabel('section Pdiff')
linkaxes(h,'x');

figure;
tiledlayout('flow');
for i = 1:length(Step)
    nexttile
    plot(PEEP,'.'); hold on;
    plot(Step(i).idx,PEEP(Step(i).idx),'r-o'); hold on; title(i);
end
end

