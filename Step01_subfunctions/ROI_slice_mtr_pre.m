function [preROIslice_filt,preBRAINslice_filt,prenuisanceslice_filt] = ROI_slice_mtr_pre(iroi,rois,islice,ROIcolumn,preROI,preBRAIN,prenuisance,preart,preexc_o,TR,order)
            fprintf('-Working on %s slice %d/%d\n',rois{iroi},islice,length(ROIcolumn));
            ROIslice = preROI(:,islice);
            BRAINslice = preBRAIN;% bc each slice has its own sets of outliers etc. can't containimate the BRAIN
            nuisanceslice = prenuisance;
            
            %remove outliers based on seed signal ( seed == each ROI slice)
            zROI = zscore(ROIslice);
            seedoutliers = find(abs(zROI)>3);
            
            ROIslice(seedoutliers,:) = NaN;
            BRAINslice(seedoutliers,:) = NaN;
            nuisanceslice(seedoutliers,:) = NaN;
            seedoutliers
            
            % remove art suspects
            if size(preexc_o,2)>6
            ROIslice(preart,:) = NaN; %censor from timeseries
            BRAINslice(preart,:) = NaN;
            nuisanceslice(preart,:) = NaN;
            end
            
            %if first or last timepoints are NaN, remove them entirely
            %(bc interpolation won't work for filtering)
            [ni nj] = find(~isnan(ROIslice));
            firstT = min(ni);
            if firstT > 1 %if there are NaNs at the beginning
                ROIslice(1:firstT-1,:) = []; %get rid of them
                BRAINslice(1:firstT-1,:) = [];
                nuisanceslice(1:firstT-1,:) = [];
            end
            [ni nj] = find(~isnan(ROIslice)); % find them again for correct indices
            lastT = max(ni);
            ROIsliceend = size(ROIslice,1);
            if lastT < ROIsliceend %if there are NaNs at the end
                ROIslice(lastT+1:end,:) = []; %get rid of them
                BRAINslice(lastT+1:end,:) = [];
                nuisanceslice(lastT+1:end,:) = [];
            end
            
            %report how many timepoints were excluded
            [iz jz] = find(isnan(ROIslice));
            fprintf('\t ---> Excluded %d timepoints, including %d at beginning and %d at end\n',...
                length(unique(iz)),firstT-1,ROIsliceend-lastT);
            
            
            %interpolate missing data for frequency filtering only
            
            %     Vq = interp1(X,V,Xq) interpolates to find Vq, the values of the
            %     underlying function V=F(X) at the query points Xq. X must
            %     be a vector of length N.
            %     If V is a vector, then it must also have length N, and Vq is the
            %     same size as Xq.  If V is an array of size [N,D1,D2,...,Dk], then
            %     the interpolation is performed for each D1-by-D2-by-...-Dk value
            %     in V(i,:,:,...,:).
            %     If Xq is a vector of length M, then Vq has size [M,D1,D2,...,Dk].
            %     If Xq is an array of size [M1,M2,...,Mj], then Vq is of size
            %     [M1,M2,...,Mj,D1,D2,...,Dk].
            
            fprintf('\t Interpolating timeseries before applying filter\n');
            ROIslice_int = ROIslice;
            ROIslice_int(isnan(ROIslice_int)) = interp1(find(~isnan(ROIslice_int)), ROIslice_int(~isnan(ROIslice_int)), find(isnan(ROIslice_int)), 'linear');
            % interpolate breaking points
            
            BRAINslice_int = BRAINslice;
            for ivox = 1:size(BRAINslice_int,2);
                tmp = BRAINslice_int(:,ivox);% work on each voxel
                tmp(isnan(tmp)) = interp1(find(~isnan(tmp)), tmp(~isnan(tmp)), find(isnan(tmp)), 'linear');% X: location of not NaN values, V:values of those non Nan
                BRAINslice_int(:,ivox) = tmp;
            end
            clear tmp
            
            nuisanceslice_int = nuisanceslice;
            for icol = 1:size(nuisanceslice_int,2);
                ntmp = nuisanceslice_int(:,icol);
                ntmp(isnan(ntmp)) = interp1(find(~isnan(ntmp)), ntmp(~isnan(ntmp)), find(isnan(ntmp)), 'linear');
                nuisanceslice_int(:,icol) = ntmp;
            end
            clear ntmp
            
            %Organized MTL & BRAINslice matrix
            ROIslice_filt = bandpass(ROIslice_int,TR,[.01 .1],order);
            BRAINslice_filt = bandpass(BRAINslice_int,TR,[.01 .1],order);
            nuisanceslice_filt = bandpass(nuisanceslice_int,TR,[.01 .1],order);%for TR=2, in mrscripts dir (change order parameter if fails)
            
            %recensor filtered data
            fprintf('\t Recensoring timeseries after applying filter\n');
            ROIslice_filt(isnan(ROIslice)) = NaN;
            BRAINslice_filt(isnan(BRAINslice)) = NaN;
            nuisanceslice_filt(isnan(nuisanceslice)) = NaN;
            

            clear BRAINslice ROIslice nuisanceslice BRAINslice_int ROIslice_int nuisanceslice_int

            %initialize vars and find non-NaN voxels
            Sp = nan(1,size(BRAINslice_filt,2));
            tmp = min(isnan(BRAINslice_filt(:,:)));
            goodvox = find(tmp==0);
         
            if isempty(goodvox)
                error('No non-NaN voxels in current timeseries.');
            end
            
            [nani nanj] = find(isnan(ROIslice_filt));
            ROIslice_filt(nani,:)=[]; %get rid of rows with NaN to speed up partialcorr
            BRAINslice_filt(nani,:)=[];
            nuisanceslice_filt(nani,:)=[];
            
            %Output
            
            preROIslice_filt = ROIslice_filt;
            preBRAINslice_filt = BRAINslice_filt;
            
            [tp v] = size(preROIslice_filt);
            [ttp vv] = size(preBRAINslice_filt);
            
            if tp == ttp
                   preruncode = [ones(tp,1) zeros(tp,1)];
            else
                error('preROIslice_filt & preBRAINslice_filt have different time points');
            end
            
            prenuisanceslice_filt = [nuisanceslice_filt preruncode];
            
end %function

            
            