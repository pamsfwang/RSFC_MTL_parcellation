    %RSFC Parcellation
    %Step02: Average correlation brains and eta map
    clear all
    close all
    clc

    disp('RSFC Parcellation- Step02: connectivity homogeneity matrices & hiearchical clustering')
    ddate = '1206';
    prepost = 'prepost';

    %Add paths
    addpath('/Users/sfwang/Desktop/Smooth_Maureen_data/Data/shared_functions/');
    addpath('/Volumes/Pam HD/Applications/spm8/');

    %Dierctory information
    outputDir = '/Users/sfwang/Desktop/Smooth_Maureen_data/Data/';
    resultsDir = ['/Volumes/Connectome_parcellation/BICRS_linkage_etadis_permutation_tests/'];
    datadirr = [outputDir,'s'];%subject directory
    maskdatadir = [outputDir,'masks/'];%masks directory: ROI masks, whole-masks ,etc

    predirectory1 = '/rest_epi_0004/';
    predirectory2= '/rest_epi_0003/';

    %Subjects
    isub = [1:6];
    numsub = length(isub);

    %File names
    rois = {'rLHC' 'rRHC' 'rLPRC_PHC' 'rRPRC_PHC'};%names of the ROI masks
    wholebrain = 'wc1_0-7_average.nii';
    %Figure color scale
    clim = [ 0.45 1];

    %% Start
    % Loading Gray matter mask
    disp('Loading whole brain mask')
    cd(maskdatadir)
    gray_header = spm_vol(wholebrain);
    gray = spm_read_vols(gray_header);
    clear gray_header

    display('Loading Correlation Maps')
    for sub = isub
        fprintf('Subjects %d \n', sub)
        %Data dirr
        dirr = [datadirr,num2str(sub,'%02d'),predirectory2];
        cd(dirr)

        %Loop through ROIs
        for iroi = 1:length(rois)

            fprintf('Loading %s \n', rois{iroi})

            ROIname = [rois{iroi},'*.nii'];
            list = dir([dirr,ROIname]);
            [rl cl] = size(list);

            for el = 1:10 % work on each slice
                name = [rois{iroi},'_s',num2str(sub,'%02d'),'_slice',num2str(el,'%02d'),'_graymatter_FC10_float32_',ddate,prepost,'.nii'];
                headerinf= spm_vol(name);
                mtr = spm_read_vols(headerinf);
                if iroi==1
                    L(sub).HPCslicemtr(:,:,:,el)= mtr;
                elseif iroi ==2
                    R(sub).HPCslicemtr(:,:,:,el)=mtr;
                elseif iroi ==3
                    L(sub).LPPslicemtr(:,:,:,el)= mtr;
                else 
                    R(sub).RPPslicemtr(:,:,:,el)=mtr;
             
                end

                numslice(iroi,sub)=el;

                clear name headerinf mtr
            end

            clear ROIname list rl
        end %iroi
    end
    disp('Finished loading correlation maps')

    disp('Check slice numbers')
    numslice

    disp('Computing average correlation brain')
    LHPCslice = numslice(1,1);
    RHPCslice = numslice(2,1);

    LPPslice = numslice(3,1);
    RPPslice = numslice(4,1);

    % Collecting correlation brains across subjects for each ROI slice
    for sub = isub;
        % Hippocampus
        for ll = 1:LHPCslice;
            tl(ll).L_HPCslice_brains(:,:,:,sub) = L(sub).HPCslicemtr(:,:,:,ll);
        end
        
        for lr = 1:RHPCslice;
            tr(lr).R_HPCslice_brains(:,:,:,sub) = R(sub).HPCslicemtr(:,:,:,lr);
        end
        
        % PRC_PHC
        for lpp = 1:LPPslice
            tlpp(lpp).L_PRCPHC_slice_brains(:,:,:,sub)= L(sub).LPPslicemtr(:,:,:,lpp);
        end
        
        for rpp = 1:RPPslice
            trpp(rpp).R_PRCPHC_slice_brains(:,:,:,sub)= R(sub).RPPslicemtr(:,:,:,rpp);
        end
        
    end

        %% Averaging correlation brains across subjects for each slice
        % Hippocmapus
        fprintf('Average correlation maps across subjects. Total number of subjects %d',numsub)
        for ll = 1:LHPCslice;
            sl(ll).slicesum = nansum(tl(ll).L_HPCslice_brains,4);
            sl(ll).sliceavg = sl(ll).slicesum/numsub;
            sl(ll).sliceavg(gray==0)=[];
            [p1 p2 p3] = size(sl(1).sliceavg);
            LHPC_avg_slice_correlation(:,ll) = sl(ll).sliceavg;%reshape(sl(ll).sliceavg,1,p1*p2*p3);
        end

        for lr = 1:RHPCslice;
            sr(lr).slicesum = nansum(tr(lr).R_HPCslice_brains,4);
            sr(lr).sliceavg = sr(lr).slicesum/numsub;
            sr(lr).sliceavg(gray==0)=[];
            [pp1 pp2 pp3] = size(sr(1).sliceavg);
            RHPC_avg_slice_correlation(:,lr) = sr(lr).sliceavg;%reshape(sr(lr).sliceavg,1,pp1*pp2*pp3);
        end

        % PHG
        for lpp= 1:LPPslice
            slpp(lpp).slicesum = nansum(tlpp(lpp).L_PRCPHC_slice_brains,4);
            slpp(lpp).sliceavg = slpp(lpp).slicesum/numsub;
            slpp(lpp).sliceavg(gray==0)=[];
            [r1 r2 r3] = size(slpp(1).sliceavg);
            LPP_avg_slice_correlation(:,lpp) = slpp(lpp).sliceavg;%reshape(slpp(lpp).sliceavg,1,r1*r2*r3);
        end

        for rpp= 1:RPPslice
            srpp(rpp).slicesum = nansum(trpp(rpp).R_PRCPHC_slice_brains,4);
            srpp(rpp).sliceavg = srpp(rpp).slicesum/numsub;
            srpp(rpp).sliceavg(gray==0)=[];
            [rr1 rr2 rr3]= size(srpp(1).sliceavg);
            RPP_avg_slice_correlation(:,rpp) = srpp(rpp).sliceavg;%reshape(srpp(rpp).sliceavg,1,rr1*rr2*rr3);
        end

        L_avg_mtr= LHPC_avg_slice_correlation;size(L_avg_mtr)
        R_avg_mtr= RHPC_avg_slice_correlation;size(R_avg_mtr)

        LPP_avg_mtr = LPP_avg_slice_correlation;size(LPP_avg_mtr)
        RPP_avg_mtr = RPP_avg_slice_correlation;size(RPP_avg_mtr)
        
        %% Computing connectivity homogeneity
        cd (resultsDir);
        
        %Hippocampus
        HPC = [L_avg_mtr R_avg_mtr];
        [hpcr hpcc] = size(HPC);

        HClnk = linkage(HPC','average','correlation');%UPGMA clustering algorithm. distance: correlation coefficients
        HCname = fullfile(resultsDir,['LR_HC_lnk']);
        %save(HCname,'HClnk')
        figure(1);dendrogram(HClnk)
        
        
        HCdis=pdist(HPC','correlation');%(1-r)
        hcsquare = squareform(1-HCdis);%connectivity homogeneity (r)
        hcidentity = eye(size(hcsquare));
        hcsquare= hcsquare+hcidentity;
        figure(2);imagesc(hcsquare)
        clear HPC 
        %PRC_PHC
        PP = [LPP_avg_mtr RPP_avg_mtr];
        PPlnk = linkage(PP','average','correlation');
        PPname = fullfile(resultsDir,['LR_PHG_lnk']);
        %save(PPname,'PPlnk')
        figure(3);dendrogram(PPlnk)
        
        
        PPdis=pdist(PP','correlation');
        ppsquare = squareform(1-PPdis);
        ppidentity = eye(size(ppsquare));
        ppsquare= ppsquare+ppidentity;
        figure(4);imagesc(ppsquare)          
        
  
disp('-------PERMUTATION TESTS-------')
        
disp('Generating weights')
    for iii = 1:10000
        for sub = isub % number of subjects: 19
            in = randi(0:1,1); % randomly select 0 or 1

            if in ==0
                index =-1; % if it is 0, make it into -1
            else
                index= 1;
            end

            subin(1,sub)=index;
        end
        ttt(iii,:)=subin;
    end

    for iii = 1:10000
        fprintf('permutations %04d \n',iii)
        % Collecting correlation brains across subjects for each ROI slice
        for sub = isub;
            % Hippocampus
            for ll = 1:LHPCslice;
                tl(ll).L_HPCslice_brains(:,:,:,sub) = L(sub).HPCslicemtr(:,:,:,ll)*ttt(iii,sub);
            end

            for lr = 1:RHPCslice;
                tr(lr).R_HPCslice_brains(:,:,:,sub) = R(sub).HPCslicemtr(:,:,:,lr)*ttt(iii,sub);
            end

            % PRC_PHC
            for lpp = 1:LPPslice
                tlpp(lpp).L_PRCPHC_slice_brains(:,:,:,sub)= L(sub).LPPslicemtr(:,:,:,lpp)*ttt(iii,sub);
            end

            for rpp = 1:RPPslice
                trpp(rpp).R_PRCPHC_slice_brains(:,:,:,sub)= R(sub).RPPslicemtr(:,:,:,rpp)*ttt(iii,sub);
            end

        end

        %% Averaging correlation brains across subjects for each slice
        % Hippocmapus
        fprintf('Average correlation maps across subjects. Total number of subjects %d',numsub)
        for ll = 1:LHPCslice;
            sl(ll).slicesum = nansum(tl(ll).L_HPCslice_brains,4);
            sl(ll).sliceavg = sl(ll).slicesum/numsub;
            sl(ll).sliceavg(gray==0)=[];
            [p1 p2 p3] = size(sl(1).sliceavg);
            LHPC_avg_slice_correlation(:,ll) = sl(ll).sliceavg;%reshape(sl(ll).sliceavg,1,p1*p2*p3);
        end

        for lr = 1:RHPCslice;
            sr(lr).slicesum = nansum(tr(lr).R_HPCslice_brains,4);
            sr(lr).sliceavg = sr(lr).slicesum/numsub;
            sr(lr).sliceavg(gray==0)=[];
            [pp1 pp2 pp3] = size(sr(1).sliceavg);
            RHPC_avg_slice_correlation(:,lr) = sr(lr).sliceavg;%reshape(sr(lr).sliceavg,1,pp1*pp2*pp3);
        end

        % PRC_PHC
        for lpp= 1:LPPslice
            slpp(lpp).slicesum = nansum(tlpp(lpp).L_PRCPHC_slice_brains,4);
            slpp(lpp).sliceavg = slpp(lpp).slicesum/numsub;
            slpp(lpp).sliceavg(gray==0)=[];
            [r1 r2 r3] = size(slpp(1).sliceavg);
            LPP_avg_slice_correlation(:,lpp) = slpp(lpp).sliceavg;%reshape(slpp(lpp).sliceavg,1,r1*r2*r3);
        end

        for rpp= 1:RPPslice
            srpp(rpp).slicesum = nansum(trpp(rpp).R_PRCPHC_slice_brains,4);
            srpp(rpp).sliceavg = srpp(rpp).slicesum/numsub;
            srpp(rpp).sliceavg(gray==0)=[];
            [rr1 rr2 rr3]= size(srpp(1).sliceavg);
            RPP_avg_slice_correlation(:,rpp) = srpp(rpp).sliceavg;%reshape(srpp(rpp).sliceavg,1,rr1*rr2*rr3);
        end

        L_avg_mtr= LHPC_avg_slice_correlation;size(L_avg_mtr)
        R_avg_mtr= RHPC_avg_slice_correlation;size(R_avg_mtr)

        LPP_avg_mtr = LPP_avg_slice_correlation;size(LPP_avg_mtr)
        RPP_avg_mtr = RPP_avg_slice_correlation;size(RPP_avg_mtr)

        %% Computing distance
        %cd (resultsDir);
        %Hippocampus
        HPC = [L_avg_mtr R_avg_mtr];
        [hpcr hpcc] = size(HPC);

        HClnkper = linkage(HPC','average','correlation');
        HCname = fullfile(resultsDir,['LR_HC_lnk_permutation_tests_',num2str(iii,'%04d')]);
        %save(HCname,'HClnkper')
        
        HCdisper = pdist(HPC','correlation');
        HCmeandis(1,iii)=mean(HCdisper);%construct null distribution
        hcsquare = squareform(1-HCdisper);
        hcidentity = eye(size(hcsquare));
        hcsquare= hcsquare+hcidentity;
        figure(100);imagesc(hcsquare)
        clear HPC HClnkper

        %PRC_PHC
        PP = [LPP_avg_mtr RPP_avg_mtr];
        PPlnkper = linkage(PP','average','correlation');
        PPname = fullfile(resultsDir,['LR_PHG_lnk_permutation_tests',num2str(iii,'%04d')]);
        %save(PPname,'PPlnkper')
        
        PPdisper = pdist(PP','correlation');
        PPmeandis(1,iii)=mean(PPdisper); %construct null distribution
        ppsquare = squareform(1-PPdisper);
        ppidentity = eye(size(ppsquare));
        ppsquare= ppsquare+ppidentity;
        figure(200);imagesc(ppsquare) 
        
        clear PP PPlnkper
        
        fprintf('finished %04d \n',iii)
    end
    
    hcmeantitle = fullfile(resultsDir,'meanHC_dis_1000');
    save(hcmeantitle,'HCmeandis')

    prcphcmeantitle = fullfile(resultsDir,'meanPHG_dis_1000');
    save(prcphcmeantitle,'PPmeandis')

    ttile = fullfile(resultsDir,'combination_weights_subjects_1000');
    save(ttile,'ttt');

    
%% Null distribution: reverse distribution
[nn xx] = hist(HCmeandis,20); %n = y axis (frequency), x: x axis ( distance values)
bx = [0:0.05:xx(1)];
ax = [bx xx];
z = zeros(1,length(bx));
an = [z nn];
h1 = axes
bar(ax, an,'facecolor',[0.7 0.7 0.7],'edgecolor',[0.7 0.7 0.7])
set(h1,'Xdir','reverse')
title('HC permutation null distribution')
disp('HC threshold')
prctile(HCmeandis,5)
clear bx ax

[n x] = hist(PPmeandis,20); %n = y axis (frequency), x: x axis ( distance values)
bx = [0:0.05:x(1)];
ax = [bx x]
z = zeros(1,length(bx));
an = [z n]
h1 = axes
bar(ax, an,'facecolor',[0.7 0.7 0.7],'edgecolor',[0.7 0.7 0.7])
set(h1,'Xdir','reverse')
title('PHG permutation null distribution')
disp('PHG threshold')
prctile(PPmeandis,5)


    %% DENDROGRAM
    hipp_labels = {'L01' 'L02' 'L03' 'L04' 'L05' 'L06' 'L07' 'L08' 'L09' 'L10' 'R01' 'R02' 'R03' 'R04' 'R05' 'R06' 'R07' 'R08' 'R09' 'R10'};
    figure(11);hd = dendrogram(HClnk,'Orientation','left','Labels',hipp_labels,'ColorThreshold',prctile(HCmeandis,5));t = title({' L R Hippocampus'; 'UPGAM Dendrogram'})
    set(hd,'LineWidth',2) 
    set(t,'FontSize',15)% title font size
    set(gca,'FontSize',13)% nodes labels font size

    prcphc_labels = {'L01' 'L02' 'L03' 'L04' 'L05' 'L06' 'L07' 'L08' 'L09' 'L10' 'R01' 'R02' 'R03' 'R04' 'R05' 'R06' 'R07' 'R08' 'R09' 'R10'};
    figure(22);hdd = dendrogram(PPlnk,'Orientation','left','Labels',prcphc_labels,'ColorThreshold',prctile(PPmeandis,5));tt = title({' L R Parahippocampal Gyrus'; 'UPGAM Dendrogram'})
    set(hdd,'LineWidth',2)
    set(tt,'FontSize',15)% title font size
    set(gca,'FontSize',13)% nodes labels font size

    