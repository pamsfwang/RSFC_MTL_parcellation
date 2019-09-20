function [preBRAIN_o, prenoise_o,preexc_o,prewholebrainvoxels] = loading_pre_resting_data(presubdatadirr,preartfile,smooth,datadirr,subjects,WMmask,CSFmask,mask, filenamewb)
%LOADING FUNCTIONAL DATA: PRE
cd(presubdatadirr)
    fprintf('-Current data directory: %s \n',presubdatadirr)
    
    fprintf('-Loading Pre Resting State Smoothed EPI Images\n')
    p = spm_select('list', presubdatadirr, ['^',smooth,'.*\.nii']);% smooth: 3mm FWHM spatial smoothing functional files
    p = spm_vol(char(p));
    lp = length(p);
    
    for it = 1:lp%number of time series
        r(it).datamtr_timeseries = spm_read_vols(p(it));
    end
    [l,w,h]= size(r(1).datamtr_timeseries);
    totalnum = numel(r(1).datamtr_timeseries);
    fprintf('\t Size of data mtr: %d %d %d. Total elements:%d \n', l,w,h,totalnum)
    
    disp('-Creating smoothed whole brain matrix for RSFC computing')
    fprintf('\t masks: %s \n',filenamewb)
    for i = 1:lp
        r(i).data = r(i).datamtr_timeseries;
        r(i).data(mask~=1)=nan;% extract gray matter voxels
        r(i).vvector = reshape(r(i).data,1,totalnum);
        r(i).vvector = r(i).vvector(~isnan(r(i).vvector));
        preBRAIN_o(i,:)= [r(i).vvector];
    end
    prewholebrainvoxels = find(~isnan(r(1).data));% whole brain voxels
    
    % LOADING SUB SPECIFIC MASKS: WM & CSF
    disp('-Loading WM & CSF')
    noisemaskdatadir = [datadirr,num2str(subjects,'%02d'),'/mprage/'];% grey matter, white matter, CSF masks dirr
    %wm
    wm_header = spm_vol(fullfile(noisemaskdatadir,WMmask));
    wm = spm_read_vols(wm_header);
    fprintf('\t WMmask: %s \n',wm_header.fname)
    
    for i = 1:lp
        n(i).wm_ts = r(i).datamtr_timeseries;
        n(i).wm_ts(wm==0)=nan; % if the mask dimentions are different from functional data, this won't work
        n(i).wm_ts = reshape(n(i).wm_ts,1, totalnum);
        n(i).mwm_ts = nanmean(n(i).wm_ts);
        wm_t(i,1) = [n(i).mwm_ts];
    end
    
    %csf
    csf_header = spm_vol(fullfile(noisemaskdatadir,CSFmask));
    csf = spm_read_vols(csf_header);
    fprintf('\t CSFmask: %s \n',csf_header.fname)
    
    for i = 1:lp
        n(i).csf_ts = r(i).datamtr_timeseries;
        n(i).csf_ts(csf==0)=nan;
        n(i).csf_ts = reshape(n(i).csf_ts,1, totalnum);
        n(i).mcsf_ts = nanmean(n(i).csf_ts);
        csf_t(i,1) = [n(i).mcsf_ts];
    end
    
    prenoise_o = [wm_t csf_t];
    
    %GET MOTION REGRESSORS + ART SUSPECTS
    preexc_o = load(preartfile);
    
    clear wm_header csf_header wm csf wm_t csf_t n p r
end %function
    