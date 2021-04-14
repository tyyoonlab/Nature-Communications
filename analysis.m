%% load files
% clear;
% load('analysis');
% z0 = 24.45; F0 = 0.03766; A1 = 57.21; A2 = 26.14; d1 = 0.9077; d2 = 3.449; % M270 (20170913 for SNU Rm2, MJ)
% z0 = 23.85; F0 = -0.02085; A1 = 12.22; A2 = 31.96; d1 = 1.562; d2 = 2.684; % M270 (20171010 for SNU Rm2, MJ_revised)
% F0 = 0.0609; z0 = 24.1605; A1 = 410.4344; d1 = 0.4211; A2 = 46.9349; d2 = 2.5711; % M270 (20171010 for SNU Rm2, MJ_revised again)
F0 = 0.0609; z0 = 24.3; A1 = 410.4344; d1 = 0.4211; A2 = 46.9349; d2 = 2.5711; % M270 (20180608 for Rm2_after fluor setup)
% z0_tmp = 24.1105;

Rbead = 1400; T = 298;
correction_factor = 0.878; pixel_size = 80;
tmp = dir('raw*');
pths = arrayfun(@(i) {[tmp(i).name,'/']}, 1:numel(tmp)); npth = numel(pths);
% [finfo,fname,nbead,nframe,fps,Roff,ori,f,t,t2,M,R,P,F,rz,rx,ry,x,y,z,dx,dy,dz] = deal(cell(npth,1));
% nfile = zeros(npth,1);
for p = 1
    disp(p);
    finfo{p} = dir([pths{p},'r*.xls']);
    nfile(p) = numel(finfo{p});
%     [fname{p},Roff{p},ori{p},f{p},M{p},R{p},P{p},F{p},x{p},y{p},z{p},dx{p},dy{p},dz{p}] = deal(cell(nfile(p),1));
%     [nbead{p},fps{p},nframe{p}] = deal(zeros(nfile(p),1));
    for n = 1:nfile(p)
        disp([int2str(n/nfile(p)*100),'% of ',pths{p}(1:end-1),'...']);
        fname{p}{n} = finfo{p}(n).name;
        fname_motor = ['s',fname{p}{n}(2:end)];
        dat = dlmread([pths{p},fname_motor]);
        t2{p}{n} = dat(:,1);
        M{p}{n} = dat(:,2);
        F{p}{n} = F0 + A1*exp(-(z0-M{p}{n})./d1) + A2*exp(-(z0-M{p}{n})./d2);
%         F{p}{n} = F0 + A1*exp(-(z0_tmp-M{p}{n})./d1) + A2*exp(-(z0_tmp-M{p}{n})./d2);
        R{p}{n} = (dat(:,3)-floor(dat(:,3)))*360;
        P{p}{n} = dat(:,4);
        
        dat = dlmread([pths{p},fname{p}{n}]);
        nframe{p}(n) = size(dat,1);
        tmp = dlmread([pths{p},'c',fname{p}{n}(2:4),'.fps']);
        fps{p}(n) = tmp(1);
%         if size(tmp,1)>1
            Roff{p}{n} = tmp(2,:);
            ori{p}{n} = tmp(3,:);
%         end
        f{p}{n} = dat(:,1);
        dat = dat(:,2:end);
        t{p}{n} = f{p}{n}/fps{p}(n);
        nbead{p}(n) = size(dat,2)/3-1;
        % subtract xy offset
        dat(:,[1:3:end,2:3:end]) = dat(:,[1:3:end,2:3:end]) - repmat(mean(dat(31:60,[1:3:end,2:3:end]),1),[nframe{p}(n),1]);
        rx{p}{n} = dat(:,1)*pixel_size;
        ry{p}{n} = dat(:,2)*pixel_size;
        rz{p}{n} = dat(:,3);
        x{p}{n} = dat(:,4:3:end)*pixel_size;
        y{p}{n} = dat(:,5:3:end)*pixel_size;
        z{p}{n} = dat(:,6:3:end);
        dx{p}{n} = (x{p}{n}-repmat(rx{p}{n},[1,nbead{p}(n)]));
        dy{p}{n} = (y{p}{n}-repmat(ry{p}{n},[1,nbead{p}(n)]));
        dz{p}{n} = (z{p}{n}-repmat(rz{p}{n},[1,nbead{p}(n)]))*correction_factor;

        % synchronize motor data
        F{p}{n} = interp1(t2{p}{n},F{p}{n},t{p}{n}); idx = find(isnan(F{p}{n}),1,'last'); F{p}{n}(1:idx) = F{p}{n}(idx+1);
        R{p}{n} = interp1(t2{p}{n},R{p}{n},t{p}{n});
        P{p}{n} = interp1(t2{p}{n},P{p}{n},t{p}{n});
    end
end
clear dat;
save('analysis_raw1');

%% plot all results
for p = 1:npth
    maxfig(10+p); clf;
    set(gcf,'defaultlinemarkersize',1,'defaultaxesfontsize',15);
    nrow = floor(sqrt(nfile(p))); ncol = ceil(nfile(p)/nrow);
    for n = 1:nfile(p)
        subplot(nrow,ncol,n);
        plot(dz{p}{n},F{p}{n},'-','linewidth',.1);
        xlim([-300, 200]); ylim([0,21]);
        title(['#',num2str(n),',',regexprep(fname{p}{n},'.xls','')]);
        grid on;
    end
    printfig(gcf,['all results_',pths{p}(1:end-1)],'png');
end
