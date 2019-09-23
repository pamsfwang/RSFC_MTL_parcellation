# RSFC_MTL_parcellation
Functional connectivity based parcellation of the human medial temporal lobe  
Wang et al., 2016

Regional differences in large-scale connectivity have been proposed to underlie functional specialization along the anterior–posterior axis of the medial temporal lobe (MTL), including the hippocampus (HC) and the parahippocampal gyrus (PHG). However, it is unknown whether functional connectivity (FC) can be used reliably to parcellate the human MTL. The current study aimed to differentiate subregions of the HC and the PHG based on patterns of whole-brain intrinsic FC. 

**Method:**   
1. FC maps were calculated for each slice along the longitudinal axis of the PHG and the HC. 
2. A hierarchical clustering algorithm was then applied to these data in order to group slices according to the similarity of their connectivity patterns.   

**Results:**  
Surprisingly, three discrete clusters were identified in the PHG. Two clusters corresponded to the parahippocampal cortex (PHC) and the perirhinal cortex (PRC), and these regions showed preferential connectivity with previously described posterior-medial and anterior-temporal networks, respectively. The third cluster corresponded to an anterior PRC region previously described as area 36d, and this region exhibited pref- erential connectivity with auditory cortical areas and with a network involved in visceral processing. The three PHG clusters showed different profiles of activation during a memory-encoding task, demonstrat- ing that the FC-based parcellation identified functionally dissociable sub-regions of the PHG. In the hip- pocampus, no sub-regions were identified via the parcellation procedure. These results indicate that connectivity-based methods can be used to parcellate functional regions within the MTL, and they suggest that studies of memory and high-level cognition need to differentiate between PHC, posterior PRC, and anterior PRC.

## Paper 
You can find the original paper of the study and the supplementary materials in this folder.

## brain_masks
In this folder, you can find three masks of bilateral MTL cortical regions (nifti format; in the MNI space) generated by the RSFC parcellation method described in the paper. Basically, the three masks are the results of the study. 
* mniPHGgroup1_PHC.nii: the first cluster identified by the RSFC parcellation method, which corresponded to the PHC
* mniPHGgroup2_postPRC.nii: the second cluster identified by the RSFC parcellation method, which corresponded to a posterior portion of the PRC
* mniPHGgroup3_antPRC.nii: the third cluster identified by the RSFC parcellation method, which corresponded to an anterior region of the PRC

Combining anterior and posterior PRC masks will give you a mask that covers the whole PRC region.

## Step01: intrinsic FC preprocessing and analysis
For each subject, functional time- series from resting-state scans were extracted from regions of interest, which were coronal slices of MTL regions in the study, and from a group-averaged mask of gray matter, individual masks of white matter (WM) and cerebrospinal fluid (CSF) (smoothed).     
**FC preprocessing**  
The first three scans in each resting-state scan were removed to allow for T1 equilibration effects. Time-series were corrected for linear trends. Time-points that were seed outliers or suspected of motion contamination were scrubbed from the time series. Seed outliers were time-points in which the seed (each coronal slice) signal deviated more than three SD from its mean. Other suspect time-points were identified via the artifact detection tools (http://www.nitrc.org/projects/arti- fact_detect), defined as time-points marked by greater than 1 mm in movement or 2% global mean signal change. Data at these time-points were interpolated and then the time-series were band pass filtered for frequencies of .01 to .1 Hz. After band pass filter- ing, the interpolated time-points were removed from the time series.     
**FC calculation**  
For each ROI coronal slice (i.e. the ROI in this study), pair-wise correlations (Pearson’s r) were computed to correlate mean time-series of each coronal slice with the time-series of all the gray matter voxels over the entire brain controlling for WM mean time-series, CSF mean time-series, six motion parameters, and session means.  

The study included two scans: "pre" indicates the first resting-state scan and "post" indicates the second resting-state scan.

## Step02: hierarchical clustering and dendrograms
**Hierarchical clustering**  
First, group-level connectivity maps were generated by averaing the FC map for each slice across subjects. Next, We used the hierarchical clustering algorithm, UPGMA (Unweighted Pair Group Method with Arithmetic Mean), implemented in MATLAB, to successively merge clusters of the coronal slices based on similarities among their FC maps. Connectivity distance was calculated for each coronal slice pair by one minus connectivity similarity (1 - r). Connectivity distance was entered into the hierarchical clustering algorithm.     

**Permutation test**  
To determine significant clusters in the dendrogram (results from the hierarchical clustering), a connectivity distance threshold was calculated via permutation tests. The null hypothesis was that the FC patterns for the coronal slices in the HC or PHG (i.e. ROIs in this study) were not heterogeneous enough to separate the coronal slices into different functional regions, and thus, there were no sub-clusters in the HC or the PHG. Permutation tests were used to determine a connectivity distance threshold at which the dendrogram was partitioned into disjoint clusters. If the connectivity distance between two sets of slices was under this threshold, this would mean that their FC patterns were no more dissimilar than what would be expected by chance, in this case, the two sets of slices would be grouped together into one functional cluster. In contrast, if the connectivity distance between two sets of slices was above the distance threshold, then their FC patterns were more dissimilar than what would be expected by chance, in which case the null hypothesis would be rejected and these two sets of slices would be separated into different func- tional clusters.  
