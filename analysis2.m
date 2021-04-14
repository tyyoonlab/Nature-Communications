%% FEC check
p = 1; n_list = ct2n([3,13,23,23,33,33,43,53,53,53],[1,1,1,2,1,2,1,1,2,3],fname{p}); scaler_list = [.974, 1.021, 1.002, 1.002, .992, .992, .973, 1.082, 1.082, 1.082];
t1_list = {[87, 128.4, 269.2], [16, 55.5, 161.9], [11, 48.8, 206.2], [126, 164.7, 264.5], [4.1, 42.8, 236.6], [8.0, 47.3, 370.7], [81.9, 120.3, 220.3], [120.7, 160.9, 262.2], [4.7, 44.1, 159.6], [74.6, 113.2, 335.5]};
t2_list = {[128.2, 167, 300], [55.5, 93, 192.5], [48.7, 86.6, 239.1], [164.7, 204, 298.6], [42.7, 82.6, 270.8], [47.3, 88.6, 405.3], [120.3, 159.2, 254.6], [160.8, 200.4, 296.3], [44.2, 83, 196.6], [113.2, 152, 367.7]};
color_list = {'SteelBlue','LightGray','OrangeRed'};

Fdat_model = .1:.1:20; T = 298;
% sophisticated handle model
Lo_PEG = 570*(5.7/80); Lp_PEG = 0.47; % bPEG (5 kD; 3.5-6.5 seems acceptable)
Lo_DNA = 1020*.338; Lp_DNA = 31; Ko_DNA = 300;
zdat_model_PEG = WLC_inv(Fdat_model,Lo_PEG,Lp_PEG,T,1);
zdat_model_DNA = eWLC_inv(Fdat_model,Lo_DNA,Lp_DNA,T,Ko_DNA,1);
% % simple handle model
% Lo_DNA = 1020*.338; Lp_DNA = 8;
% zdat_model_PEG = 0;
% zdat_model_DNA = WLC_inv(Fdat_model,Lo_DNA,Lp_DNA,T,0);
% SC model
Lp_PP = .6; nLc_PP = .4; nLc_helix = .15; % helical rise per AA
nAA_frayed_SB_HZ = 21; % 21, 24, 28 for +2, +1, 0 layer
nAA_frayed_SX_HZ = 14; % 14, 21, 28 for +4, +2, 0 layer
nAA_frayed_SB_UZ = 53; % complete
nAA_frayed_SX_UZ = 14; % 14, 21, 28 for +4, +2, 0 layer
zdat_model_FZ = zdat_model_PEG + zdat_model_DNA + 2;
zdat_model_LO = zdat_model_FZ + WLC_inv(Fdat_model,nLc_PP*(12+11),Lp_PP,T,1); % 12 and 11 from SB and SX
zdat_model_HZ = zdat_model_FZ - 2 + WLC_inv(Fdat_model,nLc_PP*(12+nAA_frayed_SB_HZ+11+nAA_frayed_SX_HZ),Lp_PP,T,1) + 2./sin(atan(2./(nAA_frayed_SB_HZ-nAA_frayed_SX_HZ)/nLc_helix));
zdat_model_UZ = zdat_model_FZ - 2 + WLC_inv(Fdat_model,nLc_PP*(12+nAA_frayed_SB_UZ+11+nAA_frayed_SX_UZ),Lp_PP,T,1) + 2./sin(atan(2./(nAA_frayed_SB_UZ-nAA_frayed_SX_UZ)/nLc_helix));
zdat_model_UF = zdat_model_FZ - 2 + WLC_inv(Fdat_model,nLc_PP*129,Lp_PP,T,1); % 65 and 64 from SB and SX

maxfig(11); clf;
set(gcf,'defaultlinelinewidth',.1);
for ni = 1:numel(n_list)
    n = n_list(ni); scaler = scaler_list(ni);
    return_dat;
    Fdat = Fdat*scaler;
    dzdat = correctFEC(Fdat,dzdat,dxdat,5,1400,Roff{p}{n},ori{p}{n});
    sel = abs(Fdat(1:end-1)-5)<.1 & diff(Fdat)>0;
    sel = sel & tdat(1:end-1)>t1_list{ni}(1) & tdat(1:end-1)<t2_list{ni}(1);
    dzdat = dzdat - mean(dzdat(sel)) + zdat_model_FZ(Fdat_model==5);
    
    subplot2(2,6,ni);
    clear h;
    for ti = 1:numel(t1_list{ni})
        [~,f1] = closest(tdat,t1_list{ni}(ti));
        [~,f2] = closest(tdat,t2_list{ni}(ti));
        frange = f1:f2; return_subdat;
        h(ti) = plot(dzdat_sub,Fdat_sub,'color',rgb(color_list{ti})); hold all;
    end
    uistack(h(1),'top')

    plot(zdat_model_FZ,Fdat_model,'-','color',rgb('blue'));
    plot(zdat_model_UZ,Fdat_model,'-','color',rgb('red'));
    plot(zdat_model_UF,Fdat_model,'-','color',rgb('gray'));
    xlim(zdat_model_FZ(Fdat_model==5)+[-50,100]);
    ylim([0,20]);
end

%%
% subplot2(2,6,11);
% p = 1; n_list = ct2n(54,1:5,fname{p}); scaler = 1.082;
% for ni = 1:numel(n_list)
%     n = n_list(ni);
%     return_dat;
%     Fdat = Fdat*scaler;
%     dzdat = correctFEC(Fdat,dzdat,dxdat,5,1400,Roff{p}{n},ori{p}{n});
%     sel = abs(Fdat(1:end-1)-5)<.1 & diff(Fdat)>0;
%     sel = sel & tdat(1:end-1)<10;
%     dzdat = dzdat - mean(dzdat(sel)) + zdat_model_FZ(Fdat_model==5);
%     dzdat = medfilt1(dzdat,12);
%     tmpf = find(Fdat>18,1,'first');
%     frange = abs(diff(Fdat(1:tmpf)))<1e-4; return_subdat;
%     plot(dzdat_sub,Fdat_sub); hold on;
% end
% plot(zdat_model_FZ,Fdat_model,'-','color',rgb('blue'));
% plot(zdat_model_UZ,Fdat_model,'-','color',rgb('red'));
% plot(zdat_model_UF,Fdat_model,'-','color',rgb('gray'));
% xlim(zdat_model_FZ(Fdat_model==5)+[-50,100]);
% ylim([0,20]);

maxfig(12); clf;
set(gcf,'defaultlinelinewidth',.1);
ni_list = [3,4,8,9,10,8,9,10];
for nii = 1:numel(ni_list)
    ni = ni_list(nii);
    if nii <= 2
        subplot(131);
    elseif nii <= 5
        subplot(132);
    else
        subplot(133);
    end
    n = n_list(ni); scaler = scaler_list(ni);
    return_dat;
    Fdat = Fdat*scaler;
    dzdat = correctFEC(Fdat,dzdat,dxdat,5,1400,Roff{p}{n},ori{p}{n});
    if nii >= 5
        dzdat = medfilt1(dzdat,5);
    end
    sel = abs(Fdat(1:end-1)-5)<.1 & diff(Fdat)>0;
    sel = sel & tdat(1:end-1)>t1_list{ni}(1) & tdat(1:end-1)<t2_list{ni}(1);
    dzdat = dzdat - mean(dzdat(sel)) + zdat_model_FZ(Fdat_model==5);
    
    clear h;
    for ti = 1
        [~,f1] = closest(tdat,t1_list{ni}(ti));
        [~,f2] = closest(tdat,t2_list{ni}(ti));
        frange = f1:f2; return_subdat;
        if ismember(ni,[3,8])
            h(ti) = plot(dzdat_sub,Fdat_sub,'color',rgb('steelblue')); hold all;
        elseif ismember(ni,[4,9])
            h(ti) = plot(dzdat_sub,Fdat_sub,'color',rgb('orangered')); hold all;
        else
            h(ti) = plot(dzdat_sub,Fdat_sub,'color',rgb('goldenrod')); hold all;
        end
    end
    uistack(h(1),'top')

    plot(zdat_model_FZ,Fdat_model,'k-');
    plot(zdat_model_LO,Fdat_model,'k-');
    plot(zdat_model_HZ,Fdat_model,'k-');
    plot(zdat_model_UZ,Fdat_model,'k-');
    plot(zdat_model_UF,Fdat_model,'k-');

    xlim([320,420]);
    ylim([5,20]);
end

%% HS FEC check
p = 1; rid_list = [63,83,83,93]; tid_list = [1,1,3,1]; scaler_list = [1.024, .976, .976, 1.014];
% t1_list = {[87, 128.4, 269.2]};
% t2_list = {[128.2, 167, 300]};
color_list = {'SteelBlue','LightGray','OrangeRed'};

Fdat_model = .1:.1:20; T = 298;
% sophisticated handle model
Lo_PEG = 570*(5.7/80); Lp_PEG = 0.47; % bPEG (5 kD; 3.5-6.5 seems acceptable)
Lo_DNA = 1020*.338; Lp_DNA = 31; Ko_DNA = 300;
zdat_model_PEG = WLC_inv(Fdat_model,Lo_PEG,Lp_PEG,T,1);
zdat_model_DNA = eWLC_inv(Fdat_model,Lo_DNA,Lp_DNA,T,Ko_DNA,1);
% % simple handle model
% Lo_DNA = 1020*.338; Lp_DNA = 8;
% zdat_model_PEG = 0;
% zdat_model_DNA = WLC_inv(Fdat_model,Lo_DNA,Lp_DNA,T,0);
% SC model
Lp_PP = .6; nLc_PP = .4; nLc_helix = .15; % helical rise per AA
nAA_frayed_SB_HZ = 21; % 21, 24, 28 for +2, +1, 0 layer
nAA_frayed_SX_HZ = 14; % 14, 21, 28 for +4, +2, 0 layer
nAA_frayed_SB_UZ = 53; % complete
nAA_frayed_SX_UZ = 14; % 14, 21, 28 for +4, +2, 0 layer
zdat_model_FZ = zdat_model_PEG + zdat_model_DNA + 2;
zdat_model_LO = zdat_model_FZ + WLC_inv(Fdat_model,nLc_PP*(12+11),Lp_PP,T,1); % 12 and 11 from SB and SX
zdat_model_HZ = zdat_model_FZ - 2 + WLC_inv(Fdat_model,nLc_PP*(12+nAA_frayed_SB_HZ+11+nAA_frayed_SX_HZ),Lp_PP,T,1) + 2./sin(atan(2./(nAA_frayed_SB_HZ-nAA_frayed_SX_HZ)/nLc_helix));
zdat_model_UZ = zdat_model_FZ - 2 + WLC_inv(Fdat_model,nLc_PP*(12+nAA_frayed_SB_UZ+11+nAA_frayed_SX_UZ),Lp_PP,T,1) + 2./sin(atan(2./(nAA_frayed_SB_UZ-nAA_frayed_SX_UZ)/nLc_helix));
zdat_model_UF = zdat_model_FZ - 2 + WLC_inv(Fdat_model,nLc_PP*129,Lp_PP,T,1); % 65 and 64 from SB and SX

maxfig(13); clf;
set(gcf,'defaultlinelinewidth',.1);
for ridi = 1:numel(rid_list)
    subplot2(2,6,ridi);
    clear h;
    for tidi = 1:2
        tid = tid_list(ridi)+(tidi-1);
        n = ct2n(rid_list(ridi),tid,fname{p}); scaler = scaler_list(ridi);
        return_dat;

        h(tidi) = plot(medfilt1(dzdat,20),Fdat); hold all;
        
        plot(zdat_model_FZ,Fdat_model,'-','color',rgb('blue'));
        plot(zdat_model_UZ,Fdat_model,'-','color',rgb('red'));
        plot(zdat_model_UF,Fdat_model,'-','color',rgb('gray'));
        xlim(zdat_model_FZ(Fdat_model==5)+[-50,100]);
        ylim([0,20]);
    end
    uistack(h(1),'top')
end

% subplot2(2,6,11);
% p = 1; n_list = ct2n(54,1:5,fname{p}); scaler = 1.082;
% for ni = 1:numel(n_list)
%     n = n_list(ni);
%     return_dat;
%     Fdat = Fdat*scaler;
%     dzdat = correctFEC(Fdat,dzdat,dxdat,5,1400,Roff{p}{n},ori{p}{n});
%     sel = abs(Fdat(1:end-1)-5)<.1 & diff(Fdat)>0;
%     sel = sel & tdat(1:end-1)<10;
%     dzdat = dzdat - mean(dzdat(sel)) + zdat_model_FZ(Fdat_model==5);
%     dzdat = medfilt1(dzdat,12);
%     tmpf = find(Fdat>18,1,'first');
%     frange = abs(diff(Fdat(1:tmpf)))<1e-4; return_subdat;
%     plot(dzdat_sub,Fdat_sub); hold on;
% end
% plot(zdat_model_FZ,Fdat_model,'-','color',rgb('blue'));
% plot(zdat_model_UZ,Fdat_model,'-','color',rgb('red'));
% plot(zdat_model_UF,Fdat_model,'-','color',rgb('gray'));
% xlim(zdat_model_FZ(Fdat_model==5)+[-50,100]);
% ylim([0,20]);

%% SCspy-a HS skewness check
p = 1; rid_list = [4,14,24,34,44,54,63,83,93]; tid_list = 1:5; scaler_list = [.974, 1.021, 1.002, .992, .973, 1.082, 1.024, .976, 1.014];
F_list = [5:14,14.5:.5:17];
cutoff_list = [76, 51, 47, 53, 58, 73, 68, 62, 50];

% F_list = setdiff(F_list,10);
dzdat_all = cell(numel(rid_list),numel(F_list));
skewval = nan(numel(rid_list),numel(F_list));
skewval_ind = nan(numel(rid_list),numel(F_list),numel(tid_list));
for ridi = 1:numel(rid_list)
% for ridi = 1
    scaler = scaler_list(ridi);
    cutoff = cutoff_list(ridi);
%     maxfig(20+ridi); clf;
%     set(gcf,'defaultlinelinewidth',.1);
    if ridi >= 7
        tid_list = 11:15;
    end
    for tidi = 1:5
        n = ct2n(rid_list(ridi),tid_list(tidi),fname{p});
        if ~isnan(n)
            return_dat;
            Fdat = Fdat*scaler;
            [dzdat,offset] = correctFEC(Fdat,dzdat,dxdat,5,1400,Roff{p}{n},ori{p}{n});
            dzdat = dzdat - offset;
            dzdat = medfilt1(dzdat,11);
%             subplot2(numel(tid_list),1,tidi);
            for Fi = 1:numel(F_list)
                [~,~,frange] = extractTrace(Fdat,dzdat,F_list(Fi),.1); return_subdat;
                frange(dzdat_sub>cutoff | abs(tdat_sub-median(tdat_sub)) > 3*std(tdat_sub) | abs(dzdat_sub)>1e3) = [];
                if F_list(Fi) == 5 || F_list(Fi) == 10
                    frange = frange(1:1200*3);
                end
                return_subdat;
                dzdat_all{ridi,Fi} = [dzdat_all{ridi,Fi}; dzdat_sub];
                skewval_ind(ridi,Fi,tidi) = skewness(dzdat_sub);
%                 plot(tdat_sub,dzdat_sub); hold all;
            end
        end
    end
    
    for Fi = 1:numel(F_list)
        dzdat_filt = dzdat_all{ridi,Fi};
%         tmp = (tmp-mean(tmp))/std(tmp);
%         [h,pval(ridi,Fi)] = kstest(tmp);
%         disp(['F = ',num2str(F_list(Fi)),' pN: ',num2str(pval(ridi,Fi))]);
        skewval(ridi,Fi) = skewness(dzdat_filt);
    end
end

% maxfig(24); clf;
% dzbin = -20:.5:100;
% Fi_list = 5:10;
% for ridi = 1:numel(rid_list)
%     for Fii = 1:numel(Fi_list)
%         subplot2(numel(rid_list),numel(Fi_list),numel(Fi_list)*(ridi-1)+Fii);
%     	histogram(dzdat_all{ridi,Fi_list(Fii)},dzbin,'linestyle','none','norm','prob');
%         xlim([0,60]);
%         ylim([0,.08]);
%         axis off;
%         if ridi == 1
%             title(F_list(Fi_list(Fii)));
%         end
%     end
% end

figure(25); clf;
for ridi = 1:numel(rid_list)
    subplot2(2,6,ridi);
%     semilogy(F_list,pval(ridi,:));	
%     plot(F_list,skewval(ridi,:),'o-'); hold all;
    ydat = nanmean(skewval_ind(ridi,:,:),3);
    eydat = nanstd(skewval_ind(ridi,:,:),[],3)/sqrt(numel(tid_list));
    errorbar(F_list,ydat,eydat,'x-'); hold all;
    hline(0);
    xlim([8,17]);
    ylim([-2,2]);
%     axis off;
    title(rid_list(ridi));
end
% legend(num2str(rid_list'),'location','nw');

%% compare HS trace before and after aSNAP
figure(13); clf;
% p = 1; n_list = ct2n(54,1:10,fname{p}); scaler = 1.082;
p = 1; n_list = ct2n(93,16:25,fname{p}); scaler = 1.014; % SC only / + aSNAP
F_list = 8:2:16;
dzdat_all = cell(2,numel(F_list));

Fdat_model = .1:.1:20; T = 298;
Lo_PEG = 570*(5/80); Lp_PEG = 0.47; % bPEG (5 kD; 3.5-6.5 seems acceptable)
Lo_dsDNA = 493*2*.338; Lp_dsDNA = 40; Ko_dsDNA = 1000;
Lo_ssDNA = 17*2*1.69*.338; Lp_ssDNA = 0.75; Ko_ssDNA = 800; % ssDNA linker
zdat_model_PEG = WLC_inv(Fdat_model,Lo_PEG,Lp_PEG,T,1);
zdat_model_dsDNA = eWLC_inv(Fdat_model,Lo_dsDNA,Lp_dsDNA,T,Ko_dsDNA,1);
zdat_model_ssDNA = mFJC(Fdat_model,Lo_ssDNA,Lp_ssDNA,T,Ko_ssDNA);

% SC model
Lp_PP = .6; nLc_PP = .4; nLc_helix = .15; % helical rise per AA
nAA_frayed_SB_HZ = 21; % 21, 24, 28 for +2, +1, 0 layer
nAA_frayed_SX_HZ = 14; % 0, 14, 21, 28 for +8, +4, +2, 0 layer
zdat_model_FZ = zdat_model_PEG + zdat_model_dsDNA + zdat_model_ssDNA + 2;
zdat_model_LO = zdat_model_FZ + WLC_inv(Fdat_model,nLc_PP*(12+11),Lp_PP,T,1); % 12 and 11 from SB and SX
zdat_model_HZ = zdat_model_FZ - 2 + WLC_inv(Fdat_model,nLc_PP*(12+nAA_frayed_SB_HZ+11+nAA_frayed_SX_HZ),Lp_PP,T,1) + 2./sin(atan(2./(nAA_frayed_SB_HZ-nAA_frayed_SX_HZ)/nLc_helix));

dzbin = 330:.5:390;
for ni = 1:numel(n_list)
    n = n_list(ni); return_dat;
    Fdat = Fdat*scaler;
    dzdat = correctFEC(Fdat,dzdat,dxdat,5,1400,Roff{p}{n},ori{p}{n});
    for Fi = 1:numel(F_list)
        [~,~,frange] = extractTrace(Fdat,dzdat,F_list(Fi),.1); return_subdat;
        if Fi == 1
            dzdat = dzdat - mean(dzdat_sub) + zdat_model_FZ(Fdat_model == F_list(Fi));
            dzdat_filt = medfilt1(dzdat,20);
        end
        return_subdat;
        dzdat_filt_sub = dzdat_filt(frange);
        
        if ni == 4
            subplot(5,6,6*(5-Fi)+(1:3));
            plot(tdat_sub-tdat_sub(1),dzdat_sub,'k','linewidth',.1); hold all;
            plot(tdat_sub-tdat_sub(1),dzdat_filt_sub,'y','linewidth',.5);
            xlim([0,1]);
            ymid = zdat_model_LO(Fdat_model==F_list(Fi));
            ylim(ymid+[-20,20]);
            hl(1) = hline(zdat_model_FZ(Fdat_model==F_list(Fi)),'b');
            hl(2) = hline(zdat_model_LO(Fdat_model==F_list(Fi)),'g');
            hl(3) = hline(zdat_model_HZ(Fdat_model==F_list(Fi)),'r'); uistack(hl(:),'bottom');
            if Fi == 1
                plot([0,0,.1]-.02,[10,0,0]+ymid-20,'k','clip','off');
                text2(-.02,ymid-20,'0.1 s','hori','left','verti','top');
                text2(-.02,ymid-20,'10 nm','hori','right','verti','bottom');
            end
%             axis off;
            text2(-.05,ymid,[num2str(F_list(Fi)),' pN'],'hori','right');
            
%             subplot(5,6,6*(5-Fi)+4);
% %             histogram(dzdat_sub,dzbin,'ori','hori','norm','prob','linestyle','none');
%             histogram(dzdat_filt_sub,dzbin,'ori','hori','norm','prob','linestyle','none');
%             xlim([0,.2]);
%             ylim(ymid+[-20,20]);
%             hl(1) = hline(zdat_model_FZ(Fdat_model==F_list(Fi)),'b');
%             hl(2) = hline(zdat_model_LO(Fdat_model==F_list(Fi)),'g');
%             hl(3) = hline(zdat_model_HZ(Fdat_model==F_list(Fi)),'r');
%             axis off;
        end
        if ni <= 5
%             dzdat_all{1,Fi} = [dzdat_all{1,Fi}; dzdat_sub];
            dzdat_all{1,Fi} = [dzdat_all{1,Fi}; dzdat_filt_sub];
        else
%             dzdat_all{2,Fi} = [dzdat_all{2,Fi}; dzdat_sub];
            dzdat_all{2,Fi} = [dzdat_all{2,Fi}; dzdat_filt_sub];
        end
    end
end

for nii = 1:2
    for Fi = 1:numel(F_list)
        subplot(5,6,6*(5-Fi)+4+(nii-1));
        dzdat_tmp = dzdat_all{nii,Fi};
        if Fi == 2 % 10 pN
            dzdat_tmp(dzdat_tmp>370) = [];
        elseif Fi == 5 % 16 pN
            dzdat_tmp(dzdat_tmp>388) = [];
        end
        histogram(dzdat_tmp,dzbin,'ori','hori','norm','prob','linestyle','none'); hold all;
        xlim([0,.2]);
        ymid = zdat_model_LO(Fdat_model==F_list(Fi));
        ylim(ymid+[-20,20]);
        hl(1) = hline(zdat_model_FZ(Fdat_model==F_list(Fi)),'b');
        hl(2) = hline(zdat_model_LO(Fdat_model==F_list(Fi)),'g');
        hl(3) = hline(zdat_model_HZ(Fdat_model==F_list(Fi)),'r');
%         axis off;
    end
end

%% disassembly data
%% extract tau values
maxfig(41); clf;
p = 1; rid_list = [93*ones(1,5), 54*ones(1,3)];
tid_list_list = {101:130, 201:230, 301:330, 331:360, 361:390, 101:150, 151:200, 201:225};
scaler_list = [1.014*ones(1,5), 1.082*ones(1,3)];
targetF_list = [16, 16, 16, 12, 8, 16, 12, 8];
dz_FZ_list = [32, 32, 32, 25, 10, 32];

% SC model
T = 298;
Lp_PP = .6; nLc_PP = .365; nLc_helix = .15; % helical rise per AA
nAA_frayed_SB_HZ = 21; % 21, 24, 28 for +2, +1, 0 layer
nAA_frayed_SX_HZ = 14; % 14, 21, 28 for +4, +2, 0 layer;  variable
nAA_frayed_SB_UZ = 53; % complete
nAA_frayed_SX_UZ = 21; % 14, 21, 28 for +4, +2, 0 layer; variable

% [taudat,tdat_all,dzdat_all] = deal(cell(numel(tid_list_list),1));
pid_list = {1:5, 6, 7, 8, 9, 10};
color_list = {'r','c','g','b','k'};
title_list = {'start','first UZ','last UZ','first UF','end'};
% for ridi = 1:5
for ridi = 6
    rid = rid_list(ridi);
    tid_list = tid_list_list{ridi};
    scaler = scaler_list(ridi);
    targetF = targetF_list(ridi);
    dz_FZ = dz_FZ_list(ridi);
    % dz_LO = dz_FZ + WLC_inv(targetF,nLc_PP*(12+11),Lp_PP,T,1); % 12 and 11 from SB and SX
    dz_HZ = dz_FZ - 2 + WLC_inv(targetF,nLc_PP*(12+nAA_frayed_SB_HZ+11+nAA_frayed_SX_HZ),Lp_PP,T,1) + 2./sin(atan(2./(nAA_frayed_SB_HZ-nAA_frayed_SX_HZ)/nLc_helix));
    dz_UZ = dz_FZ - 2 + WLC_inv(targetF,nLc_PP*(12+nAA_frayed_SB_UZ+11+nAA_frayed_SX_UZ),Lp_PP,T,1) + 2./sin(atan(2./(nAA_frayed_SB_UZ-nAA_frayed_SX_UZ)/nLc_helix));
    dz_UF = dz_FZ - 2 + WLC_inv(targetF,nLc_PP*129,Lp_PP,T,1); % 65 and 64 from SB and SX
    dz_cutoff1 = (dz_HZ+dz_UZ)/2;
    dz_cutoff2 = (dz_UZ+dz_UF)/2-1;
    if targetF < 10
        dz_cutoff2 = dz_cutoff2 - 4;
    end
    
%     taudat{ridi} = nan(numel(tid_list),3);
%     [tdat_all{ridi},dzdat_all{ridi}] = deal(nan(2401,numel(tid_list),4));
    for tidi = 7:numel(tid_list)
%     for tidi = 25
        tid = tid_list(tidi);
        n = ct2n(rid,tid,fname{p}); return_dat; Fdat = Fdat*scaler;
        if rid == 93
            if tid == 124
                dzdat = dzdat -3;
            elseif ismember(tid, [201:210])
                dzdat = dzdat +1.5;
            elseif ismember(tid, [211:230])
                dzdat = dzdat +5;
            elseif ismember(tid, [301:324])
                dzdat = dzdat +13;
            elseif ismember(tid, [325:330])
                dzdat = dzdat +12;
            elseif ismember(tid, [331:360])
                dzdat = dzdat +15;
            elseif ismember(tid, [361:376])
                dzdat = dzdat +8;
            elseif ismember(tid, [377:390])
                dzdat = dzdat +5;
            end
        elseif rid == 54
            dzdat = dzdat + 7;
        end
%             elseif ismember(tid, [305,317,323,328])
%                 dzdat = dzdat + 1;
%             end
        dzdat_filt = medfilt1(dzdat,20);

        idx = [];
        if max(Fdat) < 20 && max(Fdat) > targetF-.1
            idx = nan(4,1);
            idx(1) = find(Fdat > targetF-.1,1,'first'); % beginning of targetF
            idx(2) = find(dzdat_filt(idx(1):end) > dz_cutoff1,1,'first') + (idx(1)-1); % 1st passage to UZ
            taudat{ridi}(tidi,1) = tdat(idx(2))-tdat(idx(1)); % tau(FZ->UZ)
            dzdat_UF = dzdat_filt > dz_cutoff2;
            if rid == 93
                dzdat_UF = arrayfun(@(i) all(dzdat_UF(i+(1:600))), 0:(numel(dzdat_UF)-600));
            elseif rid == 54
                dzdat_UF = arrayfun(@(i) all(dzdat_UF(i+(1:100))), 0:(numel(dzdat_UF)-100));
            end
            tmp = find(dzdat_UF,1,'first'); % stable at UF
            if ~isempty(tmp)
                idx(4) = tmp;
                idx(3) = find(dzdat_filt(1:idx(4)) < dz_cutoff1,1,'last');
                taudat{ridi}(tidi,2) = tdat(idx(4))-tdat(idx(3)); % tau(UZ->UF)
                taudat{ridi}(tidi,3) = tdat(idx(4))-tdat(idx(1)); % tau(FZ->UF)
            end
        end
        % mark artifactual changes
        [~,idx_ignore] = findpeaks(diff(smooth(rz{p}{n},120)),'MinPeakHeight',.1,'MinPeakProminence',.1);
        dFdat = diff(smooth(Fdat,120));
        idx_ignore = [1; idx_ignore(abs(dFdat(idx_ignore))<1e-4)];

        % exclude erraneous data
        if rid == 93
            if ismember(tid,[206:210, 302,313,328,329, 378:381])
                idx = [];
            end
        end

        clf;
        if ~isempty(idx)
            pidi_end = numel(pid_list);
        else 
            pidi_end = 1;
        end
        for pidi = 1:pidi_end
            subplot(2,5,pid_list{pidi});
            yyaxis left; set(gca,'ycolor','m');
            plot(tdat,Fdat,'m');
            ylim([0,25]);
            yyaxis right; set(gca,'ycolor','k');
            plot(tdat,dzdat,'k-','linewidth',.1); hold all;
            plot(tdat,dzdat_filt,'y-','linewidth',.1);
            if pidi == 1
                xlim(tdat([1,end]));
                ylim([-40,80]);
            elseif pidi <= 5
                xlim(tdat(idx(pidi-1))+[-1,1]);
                ylim((dz_FZ-32) + [20,70]);
            else
                xlim(tdat(end)+[-2,0]);
                ylim([0,40]);
            end
            hline([dz_FZ,dz_HZ],'r-');
            hline(dz_cutoff1,'g-');
            hline(dz_cutoff2,'b-');
%             hline(dz_UF,'b--');
            hline([0,20],'k-');
            if ~isempty(idx)
                arrayfun(@(i) vline(tdat(idx(i)),[color_list{i},'-']), 1:4);
                arrayfun(@(i) text2(tdat(idx(i)),40,num2str(i),'clip','on'), 1:4);
            end
            vline(tdat(idx_ignore),'k--');
            if pidi == 1
                title(num2str([rid,tid]));
            else
                title(title_list{pidi-1},'color',color_list{pidi-1});
            end
        end
        pause;
        
        if ~isempty(idx)
            for i = 1:4
                tdat_all{ridi}(:,tidi,i) = tdat(idx(i)+(-1200:1200));
                dzdat_all{ridi}(:,tidi,i) = dzdat(idx(i)+(-1200:1200));
            end
        end
    end
end

%% compare tau values
maxfig(43); clf;
tmax_list = [12,12,12,50,250; 30,50,15,10,10; 120,120,120,120,120];
title_list = {'(-)aS,UZ','(+)aS,UZ','(+)N,UZ16','(+)N,UZ12','(+)N,UZ8';...
    '(-)aS,UF','(+)aS,UF','(+)N,UF16','(+)N,UF12','(+)N,UF8';...
    '(-)aS,tot','(+)aS,tot','(+)N,tot16','(+)N,tot12','(+)N,tot8'};
for taui = 1:3
    for ridi = 1:5
        taudat_tmp = taudat{ridi}(:,taui);
        if ridi == 3
            taudat_tmp = taudat_tmp(setdiff(1:30, [3,7,11,12,14,19,21,22,25,26,27,30]));
        end
        taudat_tmp = taudat_tmp(~(isnan(taudat_tmp) | isinf(taudat_tmp)));
        if taui == 1 || taui == 3
            taudat_tmp(taudat_tmp<=0) = [];
        end
        
        tbin = linspace(0,tmax_list(taui,ridi),8);
        subplot(3,5,5*(taui-1)+ridi);
        histogram(taudat_tmp,tbin,'norm','prob'); hold all;
        xlim(tbin([1,end]));
        ylim([0,.8]);
        set(gca,'ytick',[0:.4:.8]);

        % single exponential
        tau1 = expfit(taudat_tmp);

    %     % dobule exponential
    %     exp2pdf = @(x,p,a,b) (p)*exppdf(x,a) + (1-p)*exppdf(x,b);
    %     fitres = mle(taudat_tmp,'pdf',exp2pdf,'start',[.9,tmax_list(nii)/20,20]);

        taudat_fit = tbin+diff(tbin(1:2))/2;
%         h(1) = plot(taudat_fit,exppdf(taudat_fit,tau1)*diff(tbin(1:2))*numel(taudat_tmp),'b');
        h(1) = plot(taudat_fit,exppdf(taudat_fit,tau1)*diff(tbin(1:2)),'k');
    %     h(2) = plot(taudat_fit,exp1pdf(taudat_fit,fitres(1))*diff(tbin(1:2))*numel(taudat_tmp),'r');
    %     h(2) = plot(taudat_fit,exp2pdf(taudat_fit,fitres(1),fitres(2),fitres(3))*diff(tbin(1:2))*numel(taudat_sub),'r');
        yxlabel('# events','\tau (s)');
        title([title_list{taui,ridi},';',num2str(tau1,'%.1f'),' s']);
    %     legend(h,{'Single exp.','Double exp.'});
    end
end

%% HS intermediate
figure(43); clf;
set(gcf,'defaultlinemarkersize',5);
ridi = 4; targetF = 12;
% ridi = 3; targetF = 16;
tid_list = tid_list_list{ridi};
if ridi == 4
    rem = find(taudat{ridi}(:,1)<.1);
    rem = [rem; 10];
else
    rem = [];
end
sel_tidi = setdiff(1:numel(tid_list),rem);

% extension model
Fdat_model = .1:.1:20; T = 298;
% SC model
Lp_PP = .6; nLc_PP = .4; nLc_helix = .15; % helical rise per AA
nAA_frayed_SB_HZ = 21; % 21, 24, 28 for +2, +1, 0 layer
nAA_frayed_SX_HZ = 14; % 0, 14, 21, 28 for +8, +4, +2, 0 layer
nAA_frayed_SB_UZ = 53; % complete
nAA_frayed_SX_UZ = 21; % 14, 21, 28 for +4, +2, 0 layer
zdat_model_FZ = zeros(size(Fdat_model));
zdat_model_LO = zdat_model_FZ + WLC_inv(Fdat_model,nLc_PP*(12+11),Lp_PP,T,1); % 12 and 11 from SB and SX
zdat_model_HZ = zdat_model_FZ - 2 + WLC_inv(Fdat_model,nLc_PP*(12+nAA_frayed_SB_HZ+11+nAA_frayed_SX_HZ),Lp_PP,T,1) + 2./sin(atan(2./(nAA_frayed_SB_HZ-nAA_frayed_SX_HZ)/nLc_helix));
zdat_model_UZ = zdat_model_FZ - 2 + WLC_inv(Fdat_model,nLc_PP*(12+nAA_frayed_SB_UZ+11+nAA_frayed_SX_UZ),Lp_PP,T,1) + 2./sin(atan(2./(nAA_frayed_SB_UZ-nAA_frayed_SX_UZ)/nLc_helix));
zdat_model_UF = zdat_model_FZ - 2 + WLC_inv(Fdat_model,nLc_PP*129,Lp_PP,T,1); % 65 and 64 from SB and SX

% HMM
data = dzdat_all{ridi}(:,sel_tidi,2);
data = data - repmat(mean(data(1201+(1:60),:),1),[size(data,1),1]);
data = reshape(data,[1,size(data)]);
if ridi == 4
    mu = [-31,-28,-20,0];
else
    mu = [-45,-40,-20,0];
end
Q = numel(mu); mu = reshape(mu, [1,Q,1]);
sigma = 7*ones([1,1,Q,1]);
prior = [.8,.2,0,0];
transmat = zeros(Q);
transmat(1,2) = 1/0.1;
transmat(2,1) = 1/0.01; transmat(2,3) = 1/0.1;
transmat(3,4) = 1/0.01;
transmat = transmat*(1/1200);
for i = 1:Q
    transmat(i,i) = 1-sum(transmat(i,1:end ~= i));
end
[LL, prior, transmat, mu, sigma, mixmat] = ...
     mhmm_em(data, prior, transmat, mu, sigma, [], 'max_iter', 100, 'adj_mu', 1, 'adj_Sigma', 1, 'adj_trans', 1);

Fidx = find(Fdat_model==targetF);
tauHZ = nan(numel(sel_tidi),1);

tbin = (-100:100)*(1e-3);
dzbin = -40:2.5:10;
histmat = zeros(numel(tbin)-1,numel(dzbin)-1);
for ex = 1:numel(sel_tidi)
    tidi = sel_tidi(ex);
    xdat = tdat_all{ridi}(:,tidi,2);
    ydat = dzdat_all{ridi}(:,tidi,2); ydat = ydat - mean(ydat(1201+(1:60)));

%     % simple cutoff
%     idx1 = find(ydat(1:3600)<-25,1,'last')+1;
%     idx2 = find(ydat>-10,1,'first')-1;
    
    % HMM
    B = mixgauss_prob(data(:,:,ex), mu, sigma, mixmat);
    path = viterbi_path(prior, transmat, B);
    idx1 = find(path==3,1,'first');
    idx2 = find(path==3,1,'last');
%     
    subplot2(5,6,ex);
    sel = idx1:idx2;
    tauHZ(ex) = xdat(idx2)-xdat(idx1)+(1/1200);
    xdat = xdat - xdat(idx2+1);
    plot(xdat,ydat,'o-','linewidth',.1); hold all;
    plot(xdat(sel),ydat(sel),'ro-','linewidth',.1); hold all;
%     plot(xdat,mu(path),'k');
    xlim([-.05,.02]);
    ylim([-40,10]);

    hline(-(zdat_model_UF(Fidx)-zdat_model_UF(Fidx)),'b-');
    hline(-(zdat_model_UF(Fidx)-zdat_model_UZ(Fidx)),'b-');
    hline(-(zdat_model_UF(Fidx)-zdat_model_HZ(Fidx)),'b-');
    hline(-(zdat_model_UF(Fidx)-zdat_model_LO(Fidx)),'b-');
    hline(-(zdat_model_UF(Fidx)-zdat_model_FZ(Fidx)),'b-');
    histmat = histmat + hist2(xdat,ydat,tbin,dzbin);
end

%% location and dwell time
figure(44); clf;
subplot(121);
tspan = 10;
% imshow2(histmat,[]);
% plot(dzbin(1:end-1)+diff(dzbin(1:2))/2,histmat(101,:));
xdat = dzbin(1:end-1)+diff(dzbin(1:2))/2;
for i = 1:-1:-3
    ydat = sum(histmat(100 + i*tspan + (1:tspan),:),1);
    ydat = ydat/sum(ydat);
    plot(xdat,ydat); hold all;
end
vline(-(zdat_model_UF(Fidx)-zdat_model_UF(Fidx)),'b-');
vline(-(zdat_model_UF(Fidx)-zdat_model_UZ(Fidx)),'b-');
vline(-(zdat_model_UF(Fidx)-zdat_model_HZ(Fidx)),'b-');
vline(-(zdat_model_UF(Fidx)-zdat_model_LO(Fidx)),'b-');
vline(-(zdat_model_UF(Fidx)-zdat_model_FZ(Fidx)),'b-');
vline(mu,'k-');
yxlabel('Probability','Extension (nm)');

subplot(122);
taubin = linspace(0,25,5);
histogram(tauHZ*1e3,taubin); hold all;
tau1 = expfit(tauHZ);
taudat_fit = taubin+diff(taubin(1:2))/2;
plot(taudat_fit,exppdf(taudat_fit,tau1*1e3)*diff(taubin(1:2))*numel(tauHZ),'b');
xlim(taubin([1,end]));
% ylim([0,10]);

%% state location calculations
figure(45); clf;
% subplot(221);
% xdat = dzbin(1:end-1)+diff(dzbin(1:2))/2;
% for i = -1
%     ydat = sum(histmat(100 + i*tspan + (1:tspan),:),1);
%     ydat = ydat/sum(ydat);
%     plot(xdat,ydat); hold all;
% end
% gmfit = fittype('a1*normpdf(x,mu1,sigma1) + a2*normpdf(x,mu2,sigma2) + a3*normpdf(x,mu3,sigma3)','indep','x','coeff',{'mu1','mu2','mu3','sigma1','sigma2','sigma3','a1','a2','a3'});
% fitres = fit(xdat',ydat',gmfit,'startpoint',[-30,-20,-10,2,2,2,1,1.4,.1]);
% plot(fitres);

% extension model
T = 298; Fdat_model = 12;
Lp_PP = .6; nLc_PP = .4; nLc_helix = .15; % helical rise per AA
zdat_model_FZ = 0;
zdat_model_LO = zdat_model_FZ + WLC_inv(Fdat_model,nLc_PP*(12+11),Lp_PP,T,1); % 12 and 11 from SB and SX
zdat_model_UF = zdat_model_FZ - 2 + WLC_inv(Fdat_model,nLc_PP*129,Lp_PP,T,1); % 65 and 64 from SB and SX

% subplot(121);
% nAA_frayed_SB_HZ_list = [0,3,7,10,14,17,21,24,28]; % from +8 to 0th layer
% for i = 1:numel(nAA_frayed_SB_HZ_list)
%     nAA_frayed_SB_HZ = nAA_frayed_SB_HZ_list(i);
%     nAA_frayed_SX_HZ_list = nAA_frayed_SB_HZ_list(nAA_frayed_SB_HZ_list<=nAA_frayed_SB_HZ);
%     xdat = 8:-1:0;
%     xdat = xdat(1:numel(nAA_frayed_SX_HZ_list));
%     zdat_model_HZ = zeros(numel(nAA_frayed_SX_HZ_list),1);
%     for j = 1:numel(nAA_frayed_SX_HZ_list)
%         nAA_frayed_SX_HZ = nAA_frayed_SX_HZ_list(j);
%         zdat_model_HZ(j) = zdat_model_FZ - 2 + WLC_inv(Fdat_model,nLc_PP*(12+nAA_frayed_SB_HZ+11+nAA_frayed_SX_HZ),Lp_PP,T,1) + 2./sin(atan(2./(nAA_frayed_SB_HZ-nAA_frayed_SX_HZ)/nLc_helix));
%     end
% %     ydat = -(zdat_model_UF-zdat_model_HZ);
%     ydat = zdat_model_HZ;
%     plot(xdat,ydat,'o-'); hold all;
% end
% % hline(-(zdat_model_UF-zdat_model_LO),'g-');
% % hline(-20.8,'y-');
% hline(zdat_model_FZ,'b-');
% hline(zdat_model_LO,'g-');
% % hline((-20.8)-(-29.6529),'y-');
% set(gca,'xtick',sort(xdat),'xdir','reverse');
% yxlabel('Location (nm)','Syntaxin unfolding');

subplot(121);
nAA_frayed_SB_HZ_list = [21,24,28];
nAA_frayed_SX_HZ_list = [0,3,7,10,14,17,21];
xdat = 8:-1:2; clear ydat;
for i = 1:numel(nAA_frayed_SB_HZ_list)
    nAA_frayed_SB_HZ = nAA_frayed_SB_HZ_list(i);
    xdat = xdat(1:numel(nAA_frayed_SX_HZ_list));
    zdat_model_HZ = zeros(numel(nAA_frayed_SX_HZ_list),1);
    for j = 1:numel(nAA_frayed_SX_HZ_list)
        nAA_frayed_SX_HZ = nAA_frayed_SX_HZ_list(j);
        zdat_model_HZ(j) = zdat_model_FZ - 2 + WLC_inv(Fdat_model,nLc_PP*(12+nAA_frayed_SB_HZ+11+nAA_frayed_SX_HZ),Lp_PP,T,1) + 2./sin(atan(2./(nAA_frayed_SB_HZ-nAA_frayed_SX_HZ)/nLc_helix));
    end
%     ydat(:,i) = -(zdat_model_UF-zdat_model_HZ);    
    ydat(:,i) = zdat_model_HZ;
end
errorbar(xdat,ydat(:,2),ydat(:,2)-ydat(:,1),ydat(:,3)-ydat(:,2),'o-'); hold all;
xlim([0,8]);
ylim([0,25]);
hline(zdat_model_FZ,'b-');
hline(zdat_model_LO,'g-');
% hline((-20.8)-(-29.6529),'y-');
hline(mu(3)-mu(1),'y-');
set(gca,'xtick',sort(xdat),'xdir','reverse');

subplot(122);
% symmetric unfolding
nAA_frayed_SB_HZ_list = [0,3,7,10,14,17,21,24,28]; % from +8 to 0th layer
zdat_model_HZ = zeros(numel(nAA_frayed_SB_HZ_list),1);
for i = 1:numel(nAA_frayed_SB_HZ_list)
    nAA_frayed_SB_HZ = nAA_frayed_SB_HZ_list(i);
    nAA_frayed_SX_HZ = nAA_frayed_SB_HZ;
    zdat_model_HZ(i) = zdat_model_FZ - 2 + WLC_inv(Fdat_model,nLc_PP*(12+nAA_frayed_SB_HZ+11+nAA_frayed_SX_HZ),Lp_PP,T,1) + 2./sin(atan(2./(nAA_frayed_SB_HZ-nAA_frayed_SX_HZ)/nLc_helix));
end
xdat = 8:-1:0;
% ydat = -(zdat_model_UF-zdat_model_HZ);
ydat = zdat_model_HZ;
plot(xdat,ydat,'o-');
xlim([0,8]);
ylim([0,25]);
% hline(-(zdat_model_UF-zdat_model_LO),'g-');
% hline(-20.8,'y-');
hline(zdat_model_FZ,'b-');
hline(zdat_model_LO,'g-');
% hline((-20.8)-(-29.6529),'y-');
hline(mu(3)-mu(1),'y-');
set(gca,'xtick',sort(xdat),'xdir','reverse');

%% 16 and 12 pN disassembly example
Fdat_model = .1:.1:20; T = 298;
Lo_PEG = 570*(5/80); Lp_PEG = 0.47; % bPEG (5 kD; 3.5-6.5 seems acceptable)
Lo_dsDNA = 493*2*.338; Lp_dsDNA = 40; Ko_dsDNA = 1000;
Lo_ssDNA = 17*2*1.69*.338; Lp_ssDNA = 0.75; Ko_ssDNA = 800; % ssDNA linker
zdat_model_PEG = WLC_inv(Fdat_model,Lo_PEG,Lp_PEG,T,1);
zdat_model_dsDNA = eWLC_inv(Fdat_model,Lo_dsDNA,Lp_dsDNA,T,Ko_dsDNA,1);
zdat_model_ssDNA = mFJC(Fdat_model,Lo_ssDNA,Lp_ssDNA,T,Ko_ssDNA);
% SC model
Lp_PP = .6; nLc_PP = .4; nLc_helix = .15; % helical rise per AA
nAA_frayed_SB_HZ = 21; % 21, 24, 28 for +2, +1, 0 layer
nAA_frayed_SX_HZ = 14; % 14, 21, 28 for +4, +2, 0 layer
nAA_frayed_SB_UZ = 53; % complete
nAA_frayed_SX_UZ = 21; % 14, 21, 28 for +4, +2, 0 layer
zdat_model_FZ = zdat_model_PEG + zdat_model_dsDNA + zdat_model_ssDNA + 2;
zdat_model_LO = zdat_model_FZ + WLC_inv(Fdat_model,nLc_PP*(12+11),Lp_PP,T,1); % 12 and 11 from SB and SX
zdat_model_HZ = zdat_model_FZ - 2 + WLC_inv(Fdat_model,nLc_PP*(12+nAA_frayed_SB_HZ+11+nAA_frayed_SX_HZ),Lp_PP,T,1) + 2./sin(atan(2./(nAA_frayed_SB_HZ-nAA_frayed_SX_HZ)/nLc_helix));
zdat_model_UZ = zdat_model_FZ - 2 + WLC_inv(Fdat_model,nLc_PP*(12+nAA_frayed_SB_UZ+11+nAA_frayed_SX_UZ),Lp_PP,T,1) + 2./sin(atan(2./(nAA_frayed_SB_UZ-nAA_frayed_SX_UZ)/nLc_helix));
zdat_model_UF = zdat_model_FZ - 2 + WLC_inv(Fdat_model,nLc_PP*129,Lp_PP,T,1); % 65 and 64 from SB and SX

maxfig(91); clf;
p = 1; n_list = ct2n(93,[120,121,309,310],fname{p}); scaler = 1.014;
% p = 1; n_list = ct2n(93,[216,217,309,310],fname{p}); scaler = 1.014;
% p = 1; n_list = ct2n(93,[216,217,309,374],fname{p}); scaler = 1.014;
offset_list = [-3,-3,+8,+8];
for ni = 3:4
    subplot(2,2,ni);
    n = n_list(ni);
    return_dat;
    Fdat = Fdat*scaler;
    dzdat = correctFEC(Fdat,dzdat,dxdat,5,1400,Roff{p}{n},ori{p}{n});
    dzdat = dzdat + offset_list(ni);
    if ni == 1
        [~,f1] = closest(tdat,7.6);
        [~,f2] = closest(tdat,8.7);
        tdat = [tdat(f1-1500:f1); tdat(f2:end)-tdat(f2)+tdat(f1)]; tdat = tdat-tdat(1);
        dzdat = [dzdat(f1-1500:f1); dzdat(f2:end)];
    elseif ni == 3
        [~,f1] = closest(tdat,15.8);
        [~,f2] = closest(tdat,17);
        tdat = [tdat(f1-1500:f1); tdat(f2:end)-tdat(f2)+tdat(f1)]; tdat = tdat-tdat(1);
        dzdat = [dzdat(f1-1500:f1); dzdat(f2:end)];
    elseif ni == 2
        [~,f1] = closest(tdat,1.2);
        [~,f2] = closest(tdat,13.6);
        tdat = tdat(f1:f2)-tdat(f1);
        dzdat = dzdat(f1:f2);
    elseif ni == 4
        [~,f1] = closest(tdat,1.3);
        [~,f2] = closest(tdat,12.5);
        tdat = tdat(f1:f2)-tdat(f1);
        dzdat = dzdat(f1:f2);
    end
    dzdat_filt = medfilt1(dzdat,20);           

%         plot(tdat,dzdat,'k'); hold all;
%         plot(tdat,dzdat_filt,'y');
    nseg = ceil(numel(tdat)/18e3);
    for j = 1:nseg
        frange = (18e3*(j-1)+1):min(18e3*j,numel(tdat));
        return_subdat;
        dzdat_filt_sub = dzdat_filt(frange);
        plot(tdat_sub,dzdat_sub,'k'); hold all;
        plot(tdat_sub,dzdat_filt_sub,'y');
    end
    xlim([0,50]);
    ylim([-80,80]);

    plot([20,21],[-50,-50],'k');
    plot([20,20],[-50,-30],'k');

    if ni == 1 || ni == 3
        Fi = 160;
    else
        Fi = 20;
    end
    hline(zdat_model_UF(Fi)-zdat_model_FZ(50),'m-');
    hline(zdat_model_UZ(Fi)-zdat_model_FZ(50),'c-');
    hline(zdat_model_HZ(Fi)-zdat_model_FZ(50),'r-');
    hline(zdat_model_LO(Fi)-zdat_model_FZ(50),'g-');
    hline(zdat_model_FZ(Fi)-zdat_model_FZ(50),'b-');
    hline(zdat_model_UF(50)-zdat_model_FZ(50),'m-');
    hline(zdat_model_FZ(50)-zdat_model_FZ(50),'b-');
end

%%
maxfig(92); clf;
p = 1; n_list = ct2n(93,[337],fname{p}); scaler = 1.014;
% p = 1; n_list = ct2n(93,[216,217,309,310],fname{p}); scaler = 1.014;
% p = 1; n_list = ct2n(93,[216,217,309,374],fname{p}); scaler = 1.014;
offset_list = [9];
for ni = 1
    subplot(2,2,ni);
    n = n_list(ni);
    return_dat;
    Fdat = Fdat*scaler;
    dzdat = correctFEC(Fdat,dzdat,dxdat,5,1400,Roff{p}{n},ori{p}{n});
    dzdat = dzdat + offset_list(ni);
    if ni == 1
        [~,f1] = closest(tdat,7.8);
        [~,f2] = closest(tdat,9);
        tdat = [tdat(f1-1500:f1); tdat(f2:end)-tdat(f2)+tdat(f1)]; tdat = tdat-tdat(1);
        dzdat = [dzdat(f1-1500:f1); dzdat(f2:end)];
    end
    dzdat_filt = medfilt1(dzdat,20);           

%         plot(tdat,dzdat,'k'); hold all;
%         plot(tdat,dzdat_filt,'y');
    nseg = ceil(numel(tdat)/18e3);
    for j = 1:nseg
        frange = (18e3*(j-1)+1):min(18e3*j,numel(tdat));
        return_subdat;
        dzdat_filt_sub = dzdat_filt(frange);
        plot(tdat_sub,dzdat_sub,'k'); hold all;
        plot(tdat_sub,dzdat_filt_sub,'y');
    end
    xlim([0,50]);
    ylim([-80,80]);

    plot([20,25],[-50,-50],'k');
    plot([20,20],[-50,-30],'k');

    if ni == 1
        Fi = 120;
    end
    hline(zdat_model_UF(Fi)-zdat_model_FZ(50),'m-');
    hline(zdat_model_UZ(Fi)-zdat_model_FZ(50),'c-');
    hline(zdat_model_HZ(Fi)-zdat_model_FZ(50),'r-');
    hline(zdat_model_LO(Fi)-zdat_model_FZ(50),'g-');
    hline(zdat_model_FZ(Fi)-zdat_model_FZ(50),'b-');
    hline(zdat_model_UF(50)-zdat_model_FZ(50),'m-');
    hline(zdat_model_FZ(50)-zdat_model_FZ(50),'b-');
end

%% 20200621 Allan deviation
figure(1); clf;
p = 1; n_list = ct2n(93,11:15,fname{p}); scaler = 1.014;
% maxfig(2); clf;
% p = 1; n_list = ct2n(54,1:5,fname{p}); scaler = 1.082;
taudat = logspace(-3,0,13)';
ydat = zeros(length(taudat),numel(n_list));
ydat_filt = zeros(length(taudat),numel(n_list));
clear data;
data.rate = 1200;
for ni = 1:numel(n_list)
    n = n_list(ni);
    return_dat; Fdat = Fdat*scaler;
    [~,~,frange] = extractTrace(Fdat,dzdat,5,.1);
    frange = frange(tdat(frange) < 10);
    return_subdat;
    data.freq = dzdat_sub;
    avar = allan(data,taudat);
    ydat(:,ni) = avar.sig2;
    
    % med filtered
    data.freq = medfilt1(dzdat_sub,20);
    avar = allan(data,taudat);
    ydat_filt(:,ni) = avar.sig2;
end
errorbar(taudat,mean(ydat,2),std(ydat,[],2)/sqrt(numel(n_list))); hold all;
errorbar(taudat,mean(ydat_filt,2),std(ydat_filt,[],2)/sqrt(numel(n_list)));
set(gca,'xscale','log','yscale','log');
xlim([5e-4,2]);
ylim([.1,10]);
yxlabel('Allan deviation (nm)','Averaging time (s)');
legend({'1.2 kHz raw trace','60 Hz-filtered trace'});

%% 20200621 "Exact" location of states
Fdat_model = .1:.1:20; T = 298;
% SC model
Lp_PP = .6; nLc_PP = .4; nLc_helix = .15; % helical rise per AA
nAA_frayed_SB_HZ = 21; % 21, 24, 28 for +2, +1, 0 layer
nAA_frayed_SX_HZ = 14; % 14, 21, 28 for +4, +2, 0 layer
nAA_frayed_SB_UZ = 53; % complete
nAA_frayed_SX_UZ = 21; % 14, 21, 28 for +4, +2, 0 layer
zdat_model_FZ = zeros(size(Fdat_model));
zdat_model_LU = zdat_model_FZ + WLC_inv(Fdat_model,nLc_PP*(12+11),Lp_PP,T,1); % 12 and 11 from SB and SX
zdat_model_HU = zdat_model_FZ - 2 + WLC_inv(Fdat_model,nLc_PP*(12+nAA_frayed_SB_HZ+11+nAA_frayed_SX_HZ),Lp_PP,T,1) + 2./sin(atan(2./(nAA_frayed_SB_HZ-nAA_frayed_SX_HZ)/nLc_helix));
zdat_model_FU = zdat_model_FZ - 2 + WLC_inv(Fdat_model,nLc_PP*(12+nAA_frayed_SB_UZ+11+nAA_frayed_SX_UZ),Lp_PP,T,1) + 2./sin(atan(2./(nAA_frayed_SB_UZ-nAA_frayed_SX_UZ)/nLc_helix));
zdat_model_UC = zdat_model_FZ - 2 + WLC_inv(Fdat_model,nLc_PP*129,Lp_PP,T,1); % 65 and 64 from SB and SX

figure(2); clf;
set(gcf,'defaultlinelinewidth',5);
plot(zdat_model_FZ,Fdat_model,'color',rgb('darkblue')); hold all;
plot(zdat_model_LU,Fdat_model,'color',rgb('lightseagreen'));
plot(zdat_model_HU,Fdat_model,'color',rgb('lavender'));
plot(zdat_model_FU,Fdat_model,'color',rgb('darkseagreen'));
plot(zdat_model_UC,Fdat_model,'color',rgb('indianred'));
ylim([5,16]);
xlim([-2,37]);
grid on;
yxlabel('Force (pN)','Relative extension (nm)');

%% 20200621 aSNAP concentration
maxfig(3); clf;
p = 1; n_list = ct2n(54,[1:5,11:15,6:10],fname{p}); scaler = 1.082;
F_list = 12:2:16;
dzdat_all = cell(3,numel(F_list));

Fdat_model = .1:.1:20; T = 298;
Lo_PEG = 570*(5/80); Lp_PEG = 0.47; % bPEG (5 kD; 3.5-6.5 seems acceptable)
Lo_dsDNA = 493*2*.338; Lp_dsDNA = 40; Ko_dsDNA = 800;
Lo_ssDNA = 17*2*1.69*.338; Lp_ssDNA = 0.75; Ko_ssDNA = 800; % ssDNA linker
zdat_model_PEG = WLC_inv(Fdat_model,Lo_PEG,Lp_PEG,T,1);
zdat_model_dsDNA = eWLC_inv(Fdat_model,Lo_dsDNA,Lp_dsDNA,T,Ko_dsDNA,1);
zdat_model_ssDNA = mFJC(Fdat_model,Lo_ssDNA,Lp_ssDNA,T,Ko_ssDNA);
% SC model
Lp_PP = .6; nLc_PP = .4; nLc_helix = .15; % helical rise per AA
nAA_frayed_SB_HZ = 21; % 21, 24, 28 for +2, +1, 0 layer
nAA_frayed_SX_HZ = 14; % 0, 14, 21, 28 for +8, +4, +2, 0 layer
zdat_model_FZ = zdat_model_PEG + zdat_model_dsDNA + zdat_model_ssDNA + 2;
zdat_model_FZ = zdat_model_FZ + 14;
zdat_model_LO = zdat_model_FZ + WLC_inv(Fdat_model,nLc_PP*(12+11),Lp_PP,T,1); % 12 and 11 from SB and SX
zdat_model_HZ = zdat_model_FZ - 2 + WLC_inv(Fdat_model,nLc_PP*(12+nAA_frayed_SB_HZ+11+nAA_frayed_SX_HZ),Lp_PP,T,1) + 2./sin(atan(2./(nAA_frayed_SB_HZ-nAA_frayed_SX_HZ)/nLc_helix));

dzbin = 375:.2:405;
for ni = 1:numel(n_list)
    n = n_list(ni); return_dat;
    Fdat = Fdat*scaler;
    dzdat = correctFEC(Fdat,dzdat,dxdat,5,1400,Roff{p}{n},ori{p}{n});
    for Fi = 1:numel(F_list)
        [~,~,frange] = extractTrace(Fdat,dzdat,F_list(Fi),.1); return_subdat;
        if Fi == 1
            dzdat = dzdat - mean(dzdat_sub) + zdat_model_FZ(Fdat_model == F_list(Fi)) - 3;
            dzdat_filt = medfilt1(dzdat,20);
        end
        return_subdat;
        dzdat_filt_sub = dzdat_filt(frange);
        
        if ni <= 5
%             dzdat_all{1,Fi} = [dzdat_all{1,Fi}; dzdat_sub];
            dzdat_all{1,Fi} = [dzdat_all{1,Fi}; dzdat_filt_sub];
        elseif ni <= 10
%             dzdat_all{2,Fi} = [dzdat_all{2,Fi}; dzdat_sub];
            dzdat_all{2,Fi} = [dzdat_all{2,Fi}; dzdat_filt_sub];
        elseif ni <= 15
%             dzdat_all{2,Fi} = [dzdat_all{2,Fi}; dzdat_sub];
            dzdat_all{3,Fi} = [dzdat_all{3,Fi}; dzdat_filt_sub];
        end
    end
end

for nii = 1:3
    for Fi = 1:numel(F_list)
        subplot(3,3,3*(Fi-1)+nii);
        dzdat_tmp = dzdat_all{nii,Fi};
%         if Fi == 2 % 10 pN
%             dzdat_tmp(dzdat_tmp>370) = [];
%         elseif Fi == 5 % 16 pN
%             dzdat_tmp(dzdat_tmp>388) = [];
%         end
        histogram(dzdat_tmp,dzbin,'norm','prob','linestyle','none'); hold all;
        ylim([0,.08]);
        xlim([375,405]);
        hl(1) = vline(zdat_model_FZ(Fdat_model==F_list(Fi)),'b');
        hl(2) = vline(zdat_model_LO(Fdat_model==F_list(Fi)),'g');
        hl(3) = vline(zdat_model_HZ(Fdat_model==F_list(Fi)),'r');
%         axis off;
        title([F_list(Fi),nii]);
    end
end

%%
figure(4); clf;
n_list = [303, 314, 310];
t_offset = [.2,0,.5];
for ni = 1:numel(n_list)
    n = n_list(ni); return_dat;
    Fdat = Fdat*scaler;
    dzdat = correctFEC(Fdat,dzdat,dxdat,5,1400,Roff{p}{n},ori{p}{n});
    for Fi = 1:numel(F_list)
        [~,~,frange] = extractTrace(Fdat,dzdat,F_list(Fi),.1); return_subdat;
        if Fi == 1
            dzdat = dzdat - mean(dzdat_sub) + zdat_model_FZ(Fdat_model == F_list(Fi)) - 3;
            dzdat_filt = medfilt1(dzdat,20);
        end
        return_subdat;
        dzdat_filt_sub = dzdat_filt(frange);
    end
    subplot(3,1,ni);
    tdat_sub = tdat_sub - tdat_sub(1) - t_offset(ni);
    plot(tdat_sub, dzdat_sub,'color',rgb('lightgray')); hold all;
    plot(tdat_sub, dzdat_filt_sub,'k'); hold all;

    xlim([0,3]);
    ylim([375,405]);
    hl(1) = hline(zdat_model_FZ(Fdat_model==F_list(Fi)),'b');
    hl(2) = hline(zdat_model_LO(Fdat_model==F_list(Fi)),'g');
    hl(3) = hline(zdat_model_HZ(Fdat_model==F_list(Fi)),'r');    
end

%%
p = 1; n = ct2n(93,307,fname{p}); scaler = 1.014;
return_dat; Fdat = Fdat*scaler;

figure(201); clf;
subplot(211);
[~,f1] = closest(tdat,19.5);
[~,f2] = closest(tdat,32.5);
frange = f1:f2; return_subdat;
plot(dzdat_sub);
ylim([10,60]);


% extension model
Fdat_model = 16; T = 298;
% SC model
Lp_PP = .6; nLc_PP = .4; nLc_helix = .15; % helical rise per AA
nAA_frayed_SB_HZ = 21; % 21, 24, 28 for +2, +1, 0 layer
nAA_frayed_SX_HZ = 0; % 0, 14, 21, 28 for +8, +4, +2, 0 layer
nAA_frayed_SB_UZ = 53; % complete
nAA_frayed_SX_UZ = 14; % 14, 21, 28 for +4, +2, 0 layer
zdat_model_FZ = zeros(size(Fdat_model));
zdat_model_LO = zdat_model_FZ + WLC_inv(Fdat_model,nLc_PP*(12+11),Lp_PP,T,1); % 12 and 11 from SB and SX
zdat_model_HZ = zdat_model_FZ - 2 + WLC_inv(Fdat_model,nLc_PP*(12+nAA_frayed_SB_HZ+11+nAA_frayed_SX_HZ),Lp_PP,T,1) + 2./sin(atan(2./(nAA_frayed_SB_HZ-nAA_frayed_SX_HZ)/nLc_helix));
zdat_model_UZ = zdat_model_FZ - 2 + WLC_inv(Fdat_model,nLc_PP*(12+nAA_frayed_SB_UZ+11+nAA_frayed_SX_UZ),Lp_PP,T,1) + 2./sin(atan(2./(nAA_frayed_SB_UZ-nAA_frayed_SX_UZ)/nLc_helix));
zdat_model_UF = zdat_model_FZ - 2 + WLC_inv(Fdat_model,nLc_PP*129,Lp_PP,T,1); % 65 and 64 from SB and SX

subplot(245);
[~,f1] = closest(tdat,20.3);
[~,f2] = closest(tdat,20.8);
frange = f1:f2; return_subdat;
plot(dzdat_sub);
ylim([10,60]);
hline(zdat_model_LO+20);
hline(zdat_model_HZ+20);
hline(zdat_model_UZ+20);
hline(zdat_model_UF+20);

subplot(246);
[~,f1] = closest(tdat,31.5);
[~,f2] = closest(tdat,32);
frange = f1:f2; return_subdat;
plot(dzdat_sub);
ylim([10,60]);
hline(zdat_model_LO+20);
hline(zdat_model_HZ+20);
hline(zdat_model_UZ+20);
hline(zdat_model_UF+20);

% extension model
Fdat_model = 16; T = 298;
% SC model
% Lp_PP = .6; nLc_PP = .4; nLc_helix = .15; % helical rise per AA
Lp_PP = .77; nLc_PP = .4; nLc_helix = .15; % helical rise per AA
nAA_frayed_SB_HZ = 21; % 21, 24, 28 for +2, +1, 0 layer
nAA_frayed_SX_HZ = 0; % 0, 14, 21, 28 for +8, +4, +2, 0 layer
nAA_frayed_SB_UZ = 53; % complete
nAA_frayed_SX_UZ = 14; % 14, 21, 28 for +4, +2, 0 layer
zdat_model_FZ = zeros(size(Fdat_model));
zdat_model_LO = zdat_model_FZ + WLC_inv(Fdat_model,nLc_PP*(12+11),Lp_PP,T,1); % 12 and 11 from SB and SX
zdat_model_HZ = zdat_model_FZ - 2 + WLC_inv(Fdat_model,nLc_PP*(12+nAA_frayed_SB_HZ+11+nAA_frayed_SX_HZ),Lp_PP,T,1) + 2./sin(atan(2./(nAA_frayed_SB_HZ-nAA_frayed_SX_HZ)/nLc_helix));
zdat_model_UZ = zdat_model_FZ - 2 + WLC_inv(Fdat_model,nLc_PP*(12+nAA_frayed_SB_UZ+11+nAA_frayed_SX_UZ),Lp_PP,T,1) + 2./sin(atan(2./(nAA_frayed_SB_UZ-nAA_frayed_SX_UZ)/nLc_helix));
zdat_model_UF = zdat_model_FZ - 2 + WLC_inv(Fdat_model,nLc_PP*129,Lp_PP,T,1); % 65 and 64 from SB and SX

subplot(247);
[~,f1] = closest(tdat,20.3);
[~,f2] = closest(tdat,20.8);
frange = f1:f2; return_subdat;
plot(dzdat_sub);
ylim([10,60]);
hline(zdat_model_LO+18);
hline(zdat_model_HZ+18);
hline(zdat_model_UZ+18);
hline(zdat_model_UF+18);

subplot(248);
[~,f1] = closest(tdat,31.5);
[~,f2] = closest(tdat,32);
frange = f1:f2; return_subdat;
plot(dzdat_sub);
ylim([10,60]);
hline(zdat_model_LO+18);
hline(zdat_model_HZ+18);
hline(zdat_model_UZ+18);
hline(zdat_model_UF+18);

%% update on 20200712
p = 1; n = ct2n(93,337,fname{p}); scaler = 1.014;
return_dat; Fdat = Fdat*scaler;
[~,f0] = closest(tdat,6.7);
[~,f1] = closest(tdat,7.8);
[~,f2] = closest(tdat,8.9);
tdat(f2:end) = tdat(f2:end)-tdat(f2)+tdat(f1+1);
frange = [f0:f1,f2:nframe{p}(n)];
return_subdat;

maxfig(202); clf;
% subplot(211);
% plot(tdat_sub, dzdat_sub,'color',rgb('lightgray')); hold all;
% plot(tdat_sub, medfilt1(dzdat_sub,20),'k');

p = 1; n = ct2n(93,337,fname{p}); scaler = 1.014;
return_dat; Fdat = Fdat*scaler;
[~,f0] = closest(tdat,18);
frange = f0:nframe{p}(n);
return_subdat;
tdat_all = tdat_sub-tdat_sub(1);
dzdat_all = dzdat_sub;
Fdat_all = Fdat_sub;

p = 1; n = ct2n(93,338,fname{p}); scaler = 1.014;
return_dat; Fdat = Fdat*scaler;
[~,f0] = closest(tdat,10.2);
frange = 1:f0;
return_subdat;
tdat_all = [tdat_all; tdat_sub+tdat_all(end)+1/1200];
dzdat_all = [dzdat_all; dzdat_sub];
Fdat_all = [Fdat_all; Fdat_sub];

subplot(211);
plot(tdat_all, Fdat_all,'r'); hold all;

subplot(212);
plot(tdat_all, dzdat_all,'color',rgb('lightgray')); hold all;
plot(tdat_all, medfilt1(dzdat_all,20),'k');
