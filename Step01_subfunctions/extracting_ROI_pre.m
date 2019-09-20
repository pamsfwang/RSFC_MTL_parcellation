function [preBRAIN,preROI,prenuisance,preart, ROIcolumn] = extracting_ROI_pre(iroi,rois,roi,presubdatadirr,unsmooth,preBRAIN_o, prenoise_o,preexc_o)
    fprintf('-Working on ROI : %s \n',rois{iroi})
    cd(presubdatadirr)
    
    fprintf('-Loading Resting State unsmoothed EPI Images\n')
    np = spm_select('list', presubdatadirr, ['^',unsmooth,'.*\.nii']);%normalized functional data: wruaf
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

        preBRAIN = preBRAIN_o;
        prenoise = prenoise_o;
       
        fprintf('\t RSFC preprocessing \n')
        %DROP EARLY TIMEPOINTS
        premotion = preexc_o(4:end,1:6);
        prenoise = prenoise(4:end,:);
        
        
        preBRAIN = preBRAIN(4:end,:);
        preROI = ROI(4:end,:);
        
        %DEMEAN AND DETREND - does not work with NaN values, so do beforehand.
        premotion = spm_detrend(premotion,1);
        prenoise = spm_detrend(prenoise,1);
        
        preBRAIN = spm_detrend(preBRAIN,1);
        preROI = spm_detrend(preROI,1);
        
        %COMBINE NUISANCE VARIABLE
        prenoise = spm_orth(prenoise);
        prenuisance = [premotion prenoise];%combine motion and niose
        
        clear prenoise
        %identify art suspects
        if size(preexc_o,2)>6
            [im jm] = find(preexc_o(:,7:end)>0);
            preart = im - 3; preart = preart(preart > 0); %account for missing beginning timepoints
        preart
        else
            preart = 0;
        end
end %function
