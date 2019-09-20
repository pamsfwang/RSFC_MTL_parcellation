function [postBRAIN,postROI,postnuisance,postart, ROIcolumn] = extracting_ROI_post(iroi,rois,roi,postsubdatadirr,unsmooth,postBRAIN_o, postnoise_o,postexc_o)
    fprintf('-Working on ROI : %s \n',rois{iroi})
    cd(postsubdatadirr)
    
    fprintf('-Loading Resting State unsmoothed EPI Images\n')
    np = spm_select('list', postsubdatadirr, ['^',unsmooth,'.*\.nii']);%normalized functional data: wruaf
    np = spm_vol(char(np));
    lp = length(np);
    
    for it = 1:lp%number of time series
        nr(it).datamtr_timeseries = spm_read_vols(np(it));
    end
    [l w h] = size(nr(1).datamtr_timeseries);
        %ORGANIZE MATRIX FOR COMPUTING CORRELATIONS
        disp('-ROI time series')
        fprintf('\t ROI masks: %s\n',rois{iroi})%mask name
        for i = 1:lp
            rr(i).ROI = nr(i).datamtr_timeseries;
            rr(i).ROI(roi(iroi).ROImask==0)=nan;% roi is the structure contains ROI masks
            for yy = 1:w
                rr(i).vROI(yy,:) = reshape(rr(i).ROI(:,yy,:), 1,l*h);
            end
            rr(i).meanROI = nanmean(rr(i).vROI,2);% slice mean
            ROI(i,:)=[rr(i).meanROI];% time series data
        end
        
        
        ROIcolumn= find(~isnan(ROI(1,:)));% seed slices
        ROInancolumn = find(isnan(ROI(1,:)));% slices
        ROI(:,ROInancolumn) = [];
        fprintf('\t ROI %s slices: %d-%d total: %d \n',rois{iroi},ROIcolumn(1),ROIcolumn(end),length(ROIcolumn))
        
        clear rr nr

        postBRAIN = postBRAIN_o;
        postnoise = postnoise_o;
       
        fprintf('\t RSFC preprocessing \n')
        %DROP EARLY TIMEPOINTS
        postmotion = postexc_o(4:end,1:6);
        postnoise = postnoise(4:end,:);
        
        
        postBRAIN = postBRAIN(4:end,:);
        postROI = ROI(4:end,:);
        
        %DEMEAN AND DETREND - does not work with NaN values, so do beforehand.
        postmotion = spm_detrend(postmotion,1);
        postnoise = spm_detrend(postnoise,1);
        
        postBRAIN = spm_detrend(postBRAIN,1);
        postROI = spm_detrend(postROI,1);
        
        %COMBINE NUISANCE VARIABLE
        postnoise = spm_orth(postnoise);
        postnuisance = [postmotion postnoise];%combine motion and niose
        
        clear postnoise
        %identify art suspects
        if size(postexc_o,2)>6
            [im jm] = find(postexc_o(:,7:end)>0);
            postart = im - 3; postart = postart(postart > 0); %account for missing beginning timepoints
        postart
        else
            postart = 0;
        end
end %function
