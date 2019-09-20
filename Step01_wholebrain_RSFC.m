    %RSFC Parcellation
    %Step01: Whole-brain RSFC for each seed region
    %Author: Shao-Fang Wang, 2014

    function [] = Step01_wholebrain_RSFC()
    ddate = '1206';
    prepost = 'prepost';

    %Add paths
    addpath('/Users/sfwang/Desktop/Smooth_Maureen_data/Data/shared_functions/'); % add shared functtion directory: bandpassfilter, fisherz transformation, etc
    addpath('/Users/sfwang/Desktop/Smooth_Maureen_data/Data/Step01_subfunctions/');% add step01 sub-scripts
    addpath('/Volumes/Pam HD/Applications/spm8/');

    %Dierctory information
    outputDir = '/Users/sfwang/Desktop/Smooth_Maureen_data/Data/';
    datadirr = [outputDir,'s'];%subject directory
    maskdatadir = [outputDir,'masks/'];%masks directory: ROI masks, whole-masks ,etc

    %Subjects
    isub = [1]; % subject number

    %File names
    rois = {'rLHC' 'rRHC' 'rLPRC_PHC' 'rRPRC_PHC' 'rLPRC_PHC_overlap' 'rRPRC_PHC_overlap'};%names of the ROI masks
    filenamewb = 'wholebrainmask_rswex_18subj.nii';

    smooth = 's6wruaf'; % smoothed EPI files prefix
    unsmooth = 'wruaf'; % unsmoothed EPI files prefix
    copytempfile = 'cowruaf_float32.nii';

    preartfile = 'rp_spike_regs_pre.txt'; %artfile
    postartfile = 'rp_spike_regs_post.txt';

    WMmask = 'wc2_mask_0-8.nii';%WMmask
    CSFmask = 'wc3_mask_0-8.nii';%CSFmask
    
    %Bandpass filter
    TR = 2;
    order = 10; % 10 is default

    %% Start
    diary([outputDir,'diary_wholebrain_correlation2_',ddate,prepost]);

    %Loading masks mtr
    for iroi = 1:length(rois)
        fprintf('Loading ROI masks: %s \n',rois{iroi})
        header = spm_vol([maskdatadir,rois{iroi},'.nii']);
        roi(iroi).ROImask = spm_read_vols(header);
    end
    
    wholebrain = spm_vol(fullfile(maskdatadir,filenamewb));
    mask = spm_read_vols(wholebrain);%wholebrain mask
    
    clear header wholebrain

    for subjects = isub;
        fprintf('------------------------------------------------------------\n')
        fprintf('\n')
        fprintf('Procesing Subjects %02d \n', subjects)

        %PRE DIR
        if subjects ==10 | subjects ==14 | subjects ==17
            presubdatadirr = [datadirr,num2str(subjects,'%02d'), '/rest_epi_0004'];% resting state data folder
        else
            presubdatadirr = [datadirr,num2str(subjects,'%02d'), '/rest_epi_0003'];
        end

        %POST DIR
        if subjects ==4|subjects ==10|subjects ==14|subjects == 17;
            postsubdatadirr = [datadirr,num2str(subjects,'%02d'),'/rest_epi_0010/'];
        elseif subjects ==1 |subjects ==2;
            postsubdatadirr = [datadirr,num2str(subjects,'%02d'),'/rest_epi_0007'];
        else
            postsubdatadirr = [datadirr,num2str(subjects,'%02d'),'/rest_epi_0009/'];
        end

        fprintf('Loading Pre Data (smoothed) & Preparing first level Matrices: whole brain, noise, art file\n')
        [preBRAIN_o,prenoise_o,preexc_o,prewholebrainvoxels] = loading_pre_resting_data(presubdatadirr,preartfile,smooth,datadirr,subjects,WMmask,CSFmask,mask, filenamewb);
        fprintf('\n')
        fprintf('Loading Post Data (smoothed)& Preparing first level Matrices: whole brain, noise, art file\n')
        [postBRAIN_o, postnoise_o,postexc_o,postwholebrainvoxels] = loading_post_resting_data(postsubdatadirr,postartfile,smooth,datadirr,subjects,WMmask,CSFmask,mask, filenamewb);

        for iroi = 1:length(rois)
            fprintf('Extracting Pre ROI Data (unsmoothed) & pre-processed: whole brain, nuisance, ROI \n')
            [preBRAIN,preROI,prenuisance,preart,ROIcolumn] = extracting_ROI_pre(iroi,rois,roi,presubdatadirr,unsmooth,preBRAIN_o, prenoise_o,preexc_o);
            fprintf('\n')

            fprintf('Extracting Post ROI Data (unsmoothed) & pre-processed: whole brain, nuisance, ROI\n')
            [postBRAIN,postROI,postnuisance,postart, ROIcolumn] = extracting_ROI_post(iroi,rois,roi,postsubdatadirr,unsmooth,postBRAIN_o, postnoise_o,postexc_o);


            for islice = 1:length(ROIcolumn)% loop through each ROI slice
                fprintf('Pre ROI slice mtr: whole brai, nuisance, ROI\n')
                [preROIslice_filt,preBRAINslice_filt,prenuisanceslice_filt] = ROI_slice_mtr_pre(iroi,rois,islice,ROIcolumn,preROI,preBRAIN,prenuisance,preart,preexc_o,TR,order);
                fprintf('\n')
                fprintf('Post ROI slice mtr: whole brai, nuisance, ROI\n')
                [postROIslice_filt,postBRAINslice_filt,postnuisanceslice_filt] = ROI_slice_mtr_post(iroi,rois,islice,ROIcolumn,postROI,postBRAIN,postnuisance,postart,postexc_o,TR,order);

             
                %CALCULATE CORRELATION COEFFICIENTS
                fprintf('Concatenate PRE & POST \n')
                fprintf('-Calculating RSFC \n');
                [s1 s2] = size(preBRAINslice_filt);
                [s1post s2post] = size(postBRAINslice_filt);

                fprintf('-preBRAINslice_filt size: %d,%d; postBRAINslice_filt size: %d,%d\n',s1,s2,s1post,s2post)

                prepostROIslice_filt = [preROIslice_filt; postROIslice_filt];
                prepostBRAINslice_filt = [preBRAINslice_filt;postBRAINslice_filt];
                prepostnuisanceslice_filt = [prenuisanceslice_filt;postnuisanceslice_filt];

                [spp1 spp2]= size(prepostBRAINslice_filt);

                fprintf('-Total time points: %d; brain voxels: %d\n',spp1,spp2)
                fprintf('first and last rows of prepostnuisanceslice_filt\n')% for checking 
                prepostnuisanceslice_filt(1,:)
                prepostnuisanceslice_filt(end,:)

                for cc=1:spp2
                    cL = [prepostROIslice_filt prepostBRAINslice_filt(:,cc)];
                    SS=partialcorr(cL, prepostnuisanceslice_filt);
                    Sp(cc)=SS(1,2);
                    
                    %for tracking progress
                    if rem(cc,40000)==0
                        fprintf('%d\n',cc)
                    end
                    clear cL SS
                end


                if isequal(prewholebrainvoxels,postwholebrainvoxels)==1
                    wholebrainvoxels = prewholebrainvoxels;
                else
                    error('prewholebrainvoxels ~= postwholebrainvoxels')
                end

                fprintf('-Saving correlation map as .nii images\n')
%                 cd(presubdatadirr)
%                 name =[rois{iroi},filenamewb,'_s',num2str(subjects,'%02d'),'_slice',num2str(islice,'%02d'),'_FC10_float32_',ddate,prepost,'.nii'];% you may want to change the file name
%                 copyfile(copytempfile,name);% copy a float 32 file and use its headerinfo to save the correlation map
%                 int32fname = fullfile(presubdatadirr,name);
%                 int32headerinfo = spm_vol(int32fname);
%                 int32datamtr = spm_read_vols(int32headerinfo);
%                 int32datamtr = nan(size(int32datamtr));% clean up
                %int32datamtr(wholebrainvoxels) = Sp*10;%correlation coefficients *10 for visualization
                %spm_write_vol(int32headerinfo, int32datamtr);% save the file
                clear cL SS Sp tmp prepostBRAINslice_filt prepostROIslice_filt prepostnuisanceslice_filt name intdatamtr

                fprintf('-Completed %s slice%d\n',rois{iroi},islice)
               
                fprintf('\n')
                 
            end %ROI slice
            fprintf('Completed %s\n',rois{iroi})
             fprintf('\n')
        end %ROI
        fprintf('Completed subject: %d\n',subjects)
         fprintf('------------------------------------------------------------\n')
         fprintf('\n')

    end % subjects

    disp('Step01: whole-brian RSFC calculation completed')
    diary off
    end %function

   
