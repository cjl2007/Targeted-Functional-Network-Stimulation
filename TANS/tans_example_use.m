%% Example use of Targeted Functional Network Stimulation ("TANS")

% turn some knobs & set some options
SearchSpace = [3 18 20 24 27 28  38 53 55 59 62 63]; % these values correspond to lateral prefrontal cortex (both hemispheres)
NetworksofInterest = {'Frontoparietal'}; %
Thresholds = linspace(99.9,99,10);
CoiltoCortexDistance = 40; % in mm
GridSpacing = 2; % in mm
CoilModel = 'MagVenture_MC_B70.nii.gz';
SearchAngle = 30; % in deg.
ErrorTolerance = 5; % in mm
Sigma = 0.85; % smoothing factor, used only for visualization purposes. No need to change this. 
nThreads = 20;

% define some paths 
Paths{1} = '/home/charleslynch/SimNIBS-3.2'; % download from https://simnibs.github.io/simnibs/build/html/index.html
Paths{2} = '/home/charleslynch/MultiEchofMRI-Pipeline/res0urces/'; % this folder contains ft_read / gifti functions for reading and writing cifti files (e.g., https://github.com/MidnightScanClub/MSCcodebase). 
TansDir = pwd; % this is the dir containing all the functions & example_data folder

%% Head models

% the "tans_headmodels.m" script is a wrapper for the headreco program. In addition, a combination of FreeSurfer and Connectome workBench 
% command line utilities are used to 1) generate a skin surface and 2) calculate the distances (in geodesic space) between all pairs of skin vertices, 
% information that is used later on by other functions. 

% run the "tans_headmodels.m" module (this can take about a day to run);
tans_headmodels([TansDir '/example_data/ME01'],[TansDir '/example_data/ME01/tans'],Paths);

%% STEP 1

% the primary purpose of the "tans_roi.m" script is to identify the largest
% piece of the target network that is 1) in the search space defined by the user 
% (e.g., DLPFC) and 2) is located on a gyral crown;

% read in the functional networks for this individual;
FunctionalNetworks = ft_read_cifti_mod([TansDir '/example_data/ME01/pfm/ME01_FunctionalNetworks.dtseries.nii']);

% isolate the target network
TargetNetwork = FunctionalNetworks;
TargetNetwork.data = FunctionalNetworks.data(:,9); % note: column 9 == binary frontoparietal network map. 

% run the "tans_roi.m" module;
[TargetNetwork,LargestCluster] = tans_roi([TansDir '/example_data/ME01'],TargetNetwork,SearchSpace,[TansDir '/example_data/ME01/tans/Network_Frontoparietal'],Paths);

%% STEP 2

% the primary purpose of the "tans_searchgrid.m" script is to generate an
% array of 3D coordinates cooresponding to points on the scalp to center
% the coil center (a "search grid").

% generate a search grid on the scalp above the target network;
[SubSampledSearchGrid] = tans_searchgrid([TansDir '/example_data/ME01'],LargestCluster,[TansDir '/example_data/ME01/tans/HeadModel'],[TansDir '/example_data/ME01/tans/Network_Frontoparietal'],GridSpacing,CoiltoCortexDistance,Paths);

%% STEP 3

% the primary purpose of the "tans_simnibs.m" script is to
% run simNIBS at each point in the search grid. The results of each
% simulation are mapped to the individual's midthickness surface.

% run simNIBs at every point in the search grid;
tans_simnibs([TansDir '/example_data/ME01'],[TansDir '/example_data/ME01/tans/Network_Frontoparietal'],[TansDir '/example_data/ME01/tans/HeadModel'],SubSampledSearchGrid,CoilModel,SearchAngle,nThreads,Paths);

%% STEP 4

% run the optimization;
tans_optimize([TansDir '/example_data/ME01'],[TansDir '/example_data/ME01/tans/HeadModel'],...
[TansDir '/example_data/ME01/tans/Network_Frontoparietal'],TargetNetwork,Thresholds,Sigma,ErrorTolerance,CoilModel,Paths);

%% evaluate the extent to which the E-field "hotspot" is on target

% generate the E-field hotspot associated with the optimal
E = ft_read_cifti_mod([TansDir '/example_data/ME01/tans/Network_Frontoparietal/Optimize/normE_BestCoilCenter+BestOrientation.dtseries.nii']);
[HotSpot] = tans_make_hotspot(E,Thresholds);

% read in the surface vertex area information;
VA = ft_read_cifti_mod([TansDir '/example_data/ME01/anat/T1w/fsaverage_LR32k/ME01.midthickness_va.32k_fs_LR.dscalar.nii']);

% preallocate the on target variable;
OnTarget = zeros(length(Thresholds),...
size(FunctionalNetworks.data,2));

% sweep the E-field thresholds
for t = 1:length(Thresholds)
    
    % sweep all of the functional networks;
    for i = 1:size(FunctionalNetworks.data,2)
        ROI = FunctionalNetworks.data(1:59412,i)==1;
        OnTarget(t,i) = sum(VA.data(HotSpot(:,t+1) & ROI==1)) / sum(VA.data(HotSpot(:,t+1)==1));
    end
    
end

% load the network labels and preset color information;
NetworkLabels = dataset('XLS',[TansDir '/example_data/ME01/pfm/ME01_FunctionalNetworks.xlsx']);
NetworkColors = [NetworkLabels.R NetworkLabels.G NetworkLabels.B]/255;

H = figure; % prellocate parent figure
set(H,'position',[1 1 400 250]); hold;

% plot the horizontal stacked bar graph;
h = barh([OnTarget ; nan(1,size(OnTarget,2)) ; mean(OnTarget,1)] * 100,1,'Stacked'); hold

% change bar colors
for ii = 1:length(Networks)
    h(ii).FaceColor = [ NetworkColors(ii,1) NetworkColors(ii,2) NetworkColors(ii,3)];
end

% make it
% "pretty";
yticks(1:12);
xlim([0 100]);
ylim([0.5 12.5]);
xlabel('% of Total E-field Hotspot Surface Area');
set(gca,'FontName','Arial','FontSize',12,'TickLength',[0 0]);
yticklabels({'99.9 %','99.8 %','99.7 %','99.6 %','99.5 %','99.4 %','99.3 %','99.2 %','99.1 %','99.0 %','','Average'});

% now, to inspect the outputs, you would want to look at the following files. 

% to look at the E-field associated with the optimal coil position you can
% load the following files in wb_view

% 1) first load the wb.spec file
% e.g., wb_view ~/example_data/ME01/anat/T1w/fsaverage_LR32k/ME01.32k_fs_LR.wb.spec

% 2) then load some files that will show the outline / borders the target network
% ~/example_data/ME01/tans/Network_Frontoparietal/ROI/TargetNetwork.L.border
% ~/example_data/ME01/tans/Network_Frontoparietal/ROI/TargetNetwork.R.border

% 3) finally, load the E-field and E-field "hotspot" associatated with the optimal
% position to evaluate visually how much of the hotspot is inside the target network borders 
% ~/example_data/ME01/tans/Network_Frontoparietal/Optimize_OnTarget/normE_BestCoilCenter+BestOrientation.dtseries.nii
% ~/example_data/ME01/tans/Network_Frontoparietal/Optimize_OnTarget/normE_BestCoilCenter+BestOrientation_HotSpot.dtseries.nii

% to look at the optimal coil position on the scalp mesh, you can load the
% following. 

% 1) first load the skin mesh
% e.g., wb_view ~/example_data/ME01/tans/HeadModel/m2m_ME01/SkinSmoothed.surf.gii

% 2) then load some files that will show the coil center and handle orientation
% ~/example_data/ME01/tans/Network_Frontoparietal/Optimize_OnTarget/CoilCenter.foci
% ~/example_data/ME01/tans/Network_Frontoparietal/Optimize_OnTarget/CoilOrientation.foci
% ~/example_data/ME01/tans/Network_Frontoparietal/Optimize_OnTarget/CoilOrientation_OnTarget_s0.85.shape.gii

% say, you want to plug in the information from TANS into a neuronav program for this individual. 

% 1) first, you take the non-atlas registered T1w image as your input to the neuronav program
% ~/example_data/ME01/anat/T1w/T1w.nii.gz

% 2) the information need to set the coil center / orientation in neuronav program is contained in the following text files
% ~/example_data/ME01/tans/Network_Frontoparietal/Optimize_OnTarget/CoilCenterNative.txt
% ~/example_data/ME01/tans/Network_Frontoparietal/Optimize_OnTarget/CoilOrientationNative.txt


