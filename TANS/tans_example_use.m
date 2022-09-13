%% Example use of Targeted Functional Network Stimulation ("TANS")
% "Automated optimization of TMS coil placement for personalized functional
% network engagement" - Lynch et al., 2022 (Neuron)

% define some paths 
Paths{1} = '~/SimNIBS-3.2'; % download from https://simnibs.github.io/simnibs/build/html/index.html
Paths{2} = '~/Utilities/'; % this folder contains ft_read / gifti functions for reading and writing cifti files (e.g., https://github.com/MidnightScanClub/MSCcodebase). 
Paths{3} = '~/TANS/'; % 

% add folders 
% to search path;
for i = 1:length(Paths)
    addpath(genpath(Paths{i}));
end

% Before proceeding, confirm that FSL, FreeSurfer, Connectome Workbench,
% and SimNIBS command line utilities are available when using system commands in Matlab. 

% If successful, each of the commands below should return status as 0. 

% Confirm various software is available
[status,~] = system('mris_convert -version'); % freesurfer
[status,~] = system('headreco --version'); % simnibs
[status,~] = system('wb_command -version'); % connectome workbench
[status,~] = system('flirt -version'); % fsl

% If successful, each of the commands below should return status as 2. 

% Confirm that functions for reading and writing 
% CIFTI and GIFTI files are also available
status = exist('ft_read_cifti_mod','file'); 
status = exist('gifti','file');

% define the path to 
% directory containing example data;
DataDir = '/path/to/example_data'; % note: should be obvious, but user needs to manually define this

%% Create the head model
%
% Timing: Typically several hours to a day. 
% 
% Use the function “tans_headmodels” to create a tetrahedral head mesh for electric field 
% modeling using the headreco program (Nielsen et al., 2018). In addition, this function 
% generates skin surface geometry files (*.surf.gii) and calculates the distance in geodesic 
% space between all pairs of skin vertices. 

% The inputs to the function are:
% Subject: The subject ID (string).
% T1w: Path to T1w-weighted anatomical image (string).
% T2w: Path to T2w-weighted anatomical image (string). 
% SmoothingStrength: Amount of spatial smoothing to apply to the skin surface (possible range: 0.0 - 1.0, numeric)
% SmoothingIterations: The number of iterations when smoothing the skin surface (numeric). 
% OutDir : Path to the output folder (string).
% Paths: Paths to folders that must be added to Matlab search path (cell array of strings). 

% define inputs
Subject = 'ME01';
SmoothingStrength = 0.25;
SmoothingIterations = 500;
T1w = [DataDir '/' Subject '/anat/T1w/T1w_acpc_dc_restore.nii.gz'];
T2w = [DataDir '/' Subject '/anat/T1w/T2w_acpc_dc_restore.nii.gz'];
OutDir = [DataDir '/' Subject '/tans'];

% run the tans_headmodels function
tans_headmodels(Subject,T1w,T2w,SmoothingStrength,SmoothingIterations,OutDir,Paths);

%% Identify the target network patch. 
%
% Timing: Typically less than 1 minute.  
% 
% Functional networks are often multifocal (they tend to have representation in multiple cortical zones),
% and it can be unclear a priori which part of the network TMS should be delivered in order to maximize 
% stimulation specificity. Use the function “tans_roi” to find the largest patch of the target network 
% on a gyral crown within the search space specified by the user (for example, lateral prefrontal cortex).
% The rationale is two-fold. First, targeting the largest piece of the target network decreases the
% likelihood of inadvertently stimulating non-target brain regions or networks. Second, the E-field is 
% stronger on the gyral crown than it is in the sulcus. 

% The inputs to the function are:
% TargetNetwork: A CIFTI file containing the functional network of interest (structure array). Non-zero values in TargetNetwork.data are considered target network vertices.
% MidthickSurfs: Paths to low (32k) dimensional FS_LR midthickness surfaces (Cell array of strings, MidthickSurfs{1} = path to LH, MidthickSurfs{2} = path to RH).
% VertexSurfaceArea: A CIFTI file containing the vertex surface areas (structure array). 
% Sulc:  CIFTI file containing sulcal depth information (structure array). 
% SearchSpace: A CIFTI file containing a binarized mask representing the search space (structure array). 
% OutDir: Path to the output folder (string).
% Paths: Paths to folders that must be added to Matlab search path (cell array of strings). 

% read in the functional networks for this individual;
FunctionalNetworks = ft_read_cifti_mod([DataDir '/' Subject '/pfm/' Subject '_FunctionalNetworks.dtseries.nii']);

% isolate the target network
TargetNetwork = FunctionalNetworks;
TargetNetwork.data(TargetNetwork.data~=9) = 0; % note: 9 == frontoparietal network map. 
TargetNetwork.data(TargetNetwork.data~=0) = 1; % binarize.

% load gyral labels and create a binary search space mask;
GyralLabels = ft_read_cifti_mod([DataDir '/' Subject '/anat/MNINonLinear/fsaverage_LR32k/' Subject '.aparc.32k_fs_LR.dlabel.nii']);
SearchSpace = FunctionalNetworks; % preallocate;
SearchSpace.data = zeros(size(SearchSpace.data,1),1); % blank slate;
SearchSpace.data(1:59412,1) = ismember(GyralLabels.data,[3 18 19 20 24 27 28 38 53 54 55 59 62 63]); % note: these values correspond to bilateral lateral prefrontal cortex and were manually extracted. 

% load sulcal depth information;
Sulc = ft_read_cifti_mod([DataDir '/' Subject '/anat/MNINonLinear/fsaverage_LR32k/' Subject '.sulc.32k_fs_LR.dscalar.nii']);
BrainStructure = GyralLabels.brainstructure; % barrow the brain structure index
Sulc.data(BrainStructure==-1) = []; % remove medial wall vertices present in this file. 

% define input variables;
MidthickSurfs{1} = [DataDir '/' Subject '/anat/T1w/fsaverage_LR32k/' Subject '.L.midthickness.32k_fs_LR.surf.gii']; 
MidthickSurfs{2} = [DataDir '/' Subject '/anat/T1w/fsaverage_LR32k/' Subject '.R.midthickness.32k_fs_LR.surf.gii']; 
VertexSurfaceArea = ft_read_cifti_mod([DataDir '/' Subject '/anat/T1w/fsaverage_LR32k/' Subject '.midthickness_va.32k_fs_LR.dscalar.nii']);
OutDir = [DataDir '/' Subject '/tans/Network_Frontoparietal'];

% run the tans_roi function
tans_roi(TargetNetwork,MidthickSurfs,VertexSurfaceArea,Sulc,SearchSpace,OutDir,Paths);

% Note: The TANS approach is intended to target surface-registered functional brain networks mapped in individual subjects. 
% In principle, however, any spatial map, such as an ICA component, functional connectivity or task activation map can be used instead. 
% Weighted maps must be first thresholded and binarized before they can be used as inputs to the tans_roi function. 

% Note: We recommend visually examining the target network patch ([OutDir '/ROI/TargetNetworkPatch.dtseries.nii’]). 
% It is important to confirm that the patch falls within the search space of interest and that it appears amenable to TMS 
% before proceeding to next steps. A patch that is small (< 1000 mm2) or located on the medial or inferior surface of the
% brain is unlikely to produce good results. 

%% Make a search grid on the scalp above the target network patch centroid.
%
% Timing: Typically 5-10 minutes.
%
% Use the tans_searchgrid function to generate an array of three-dimensional coordinates 
% representing a search grid on the scalp directly above the centroid of the target network patch.
 
% The inputs to the function are:
% TargetNetworkPatch: A CIFTI file containing the target network patch (structure array). Non-zero values in TargetNetworkPatch.data are considered target network patch vertices.
% PialSurfs: Paths to low (32k) dimensional FS_LR pial surfaces (Cell array of strings, PialSurfs{1} = path to LH, PialSurfs{2} = path to RH).
% SkinVol: Path to the volumetric skin segmentation volume (string). 
% SkinSurf: Path to the skin surface geometry file (string). 
% SkinDistanceMatrix: Path to the skin distance matrix (string). This is a .m file. 
% SearchGridRadius: Radius (in mm) of the search grid on the scalp directly above the centroid of the target network patch. 
% GridSpacing: The (approximate) distance between pairs of vertices in the search grid, in mm. 
% OutDir: Path to the output folder (string).
% Paths: Paths to folders that must be added to Matlab search path (cell array of strings). 

% define inputs
Subject = 'ME01';
OutDir = [DataDir '/' Subject '/tans/Network_Frontoparietal'];
TargetNetworkPatch = ft_read_cifti_mod([OutDir '/ROI/TargetNetworkPatch.dtseries.nii']);
PialSurfs{1} = [DataDir '/' Subject '/anat/T1w/fsaverage_LR32k/' Subject '.L.pial.32k_fs_LR.surf.gii']; 
PialSurfs{2} = [DataDir '/' Subject '/anat/T1w/fsaverage_LR32k/' Subject '.R.pial.32k_fs_LR.surf.gii'];
SkinVol = [DataDir '/' Subject '/tans/HeadModel/m2m_' Subject '/skin.nii.gz']; 
SkinSurf = [DataDir '/' Subject '/tans/HeadModel/m2m_' Subject '/Skin.surf.gii'];
SkinDistanceMatrix = [DataDir '/' Subject '/tans/HeadModel/m2m_' Subject '/SkinDistanceMatrix.mat']; 
GridSpacing = 2; 
SearchGridRadius = 40; 

% run the tans_searchgrid function
[SubSampledSearchGrid,FullSearchGrid] = tans_searchgrid(TargetNetworkPatch,PialSurfs,SkinVol,SkinSurf,SkinDistanceMatrix,GridSpacing,SearchGridRadius,OutDir,Paths);

% Note: We recommend reducing the density or radius of the search grid to reduce the total number of simulations that will 
% be performed when the total run time is a factor.

%% Perform electric field modeling iteratively in the search grid.
%
% Timing: Depends on the total number of simulations performed and available parallel pools in Matlab, typically several hours.
%
% Use the tans_simnibs function to perform electric field simulations at each point in the search grid using 
% SimNibs (Thielscher et al., 2015). A number of coil orientations determined by an angle (in degrees) 
% specified by the user are performed. The strength of the electric field is mapped to the individual’s midthickness surface.

% The inputs to the function are:
% SearchGridCoords: Number of coil positions x 3 numeric array generated by tans_searchgrid function. 
% HeadMesh: Path to tetrahedral head mesh (string).
% CoilModel: Name of the coil model (string). Must be available in SimNIBS library ( /SimNIBS-3.2/simnibs/ccd-files/). 
% AngleResolution: Inter angle resolution used in the search grid in mm (numeric).
% SkinSurf: Path to the skin surface geometry file (string).
% MidthickSurfs: Paths to low (32k) dimensional FS_LR midthickness surfaces (Cell array of strings, MidthickSurfs{1} = path to LH, MidthickSurfs{2} = path to RH).
% WhiteSurfs: Paths to low (32k) dimensional FS_LR white surfaces (Cell array of strings, WhiteSurfs{1} = path to LH, WhiteSurfs{2} = path to RH).
% PialSurfs: Paths to low (32k) dimensional FS_LR pial surfaces (Cell array of strings, PialSurfs{1} = path to LH, PialSurfs{2} = path to RH).
% MedialWallMasks: Paths to low (32k) dimensional FS_LR medial wall masks (Cell array of strings, MedialWallMasks{1} = path to LH, MedialWallMasks{2} = path to RH).
% nThreads: The number of parallel pools that will be used in Matlab (numeric).
% OutDir: Path to the output folder (string).
% Paths: Paths to folders that must be added to Matlab search path (cell array of strings). 

% define inputs
SearchGridCoords = SubSampledSearchGrid;
PialSurfs{1} = [DataDir '/' Subject '/anat/T1w/fsaverage_LR32k/' Subject '.L.pial.32k_fs_LR.surf.gii']; 
PialSurfs{2} = [DataDir '/' Subject '/anat/T1w/fsaverage_LR32k/' Subject '.R.pial.32k_fs_LR.surf.gii']; 
WhiteSurfs{1} = [DataDir '/' Subject '/anat/T1w/fsaverage_LR32k/' Subject '.L.white.32k_fs_LR.surf.gii']; 
WhiteSurfs{2} = [DataDir '/' Subject '/anat/T1w/fsaverage_LR32k/' Subject '.R.white.32k_fs_LR.surf.gii']; 
MidthickSurfs{1} = [DataDir '/' Subject '/anat/T1w/fsaverage_LR32k/' Subject '.L.midthickness.32k_fs_LR.surf.gii']; 
MidthickSurfs{2} = [DataDir '/' Subject '/anat/T1w/fsaverage_LR32k/' Subject '.R.midthickness.32k_fs_LR.surf.gii']; 
MedialWallMasks{1} = [DataDir '/' Subject '/anat/MNINonLinear/fsaverage_LR32k/' Subject '.L.atlasroi.32k_fs_LR.shape.gii']; 
MedialWallMasks{2} = [DataDir '/' Subject '/anat/MNINonLinear/fsaverage_LR32k/' Subject '.R.atlasroi.32k_fs_LR.shape.gii'];
SkinSurf = [DataDir '/' Subject '/tans/HeadModel/m2m_' Subject '/Skin.surf.gii']; 
HeadMesh = [DataDir '/' Subject '/tans/HeadModel/' Subject '.msh']; 
CoilModel = 'MagVenture_MC_B70.nii.gz';  
AngleResolution = 30; 
nThreads = 20; 

% run the tans_simnibs function
tans_simnibs(SubSampledSearchGrid,HeadMesh,CoilModel,AngleResolution,SkinSurf,MidthickSurfs,WhiteSurfs,PialSurfs,MedialWallMasks,nThreads,OutDir,Paths);

% Note: A fixed stimulation intensity (𝑑𝐼/𝑑𝑡 = 1 A/µs) is used during this stage of TANS because the strength of the E-field 
% varies linearly with 𝑑𝐼/𝑑𝑡 (the speed of variation of the current throughout the coil) and, for this reason, 
% has no effect on its spatial distribution (including where it is maximal relative to the target).

%% Find the coil placement that best aligns the electric field hotspot with the target network
%
% Timing: Depends on the total number of simulations performed, typically 1-2 hours. 
%
% Use the tans_optimize function to find which of the coil placements submitted to the tans_simnibs function 
% produced the greatest on-target value (calculated as the proportion of the E-field hotspot inside the
% target network). In other words, which coil placement was best for maximizing stimulation specificity.

% The inputs to the function are:
% Subject: The subject ID (string).
% TargetNetwork: A CIFTI file containing the functional network of interest (structure array). Non-zero values in TargetNetwork.data are considered target network vertices.
% AvoidanceRegion: A CIFTI file indexing the brain regions or networks of no interest (structure array). Non-zero values in AvoidanceRegion.data are considered non-target network vertices. 
% PercentileThresholds: A range of percentiles used for operationalizing the E-field hotspot (numeric). For example, linspace(99.9,99,10).
% SkinVol: Path to the volumetric skin segmentation volume (string). 
% SkinSurf: Path to the skin surface geometry file (string). 
% SkinDistanceMatrix: Path to the skin distance matrix (string). This is a .m file. 
% VertexSurfaceArea: A CIFTI file containing the vertex surface areas (structure array). 
% MidthickSurfs: Paths to low (32k) dimensional FS_LR midthickness surfaces (Cell array of strings, MidthickSurfs{1} = path to LH, MidthickSurfs{2} = path to RH).
% WhiteSurfs: Paths to low (32k) dimensional FS_LR white surfaces (Cell array of strings, WhiteSurfs{1} = path to LH, WhiteSurfs{2} = path to RH).
% PialSurfs: Paths to low (32k) dimensional FS_LR pial surfaces (Cell array of strings, PialSurfs{1} = path to LH, PialSurfs{2} = path to RH).
% MedialWallMasks: Paths to low (32k) dimensional FS_LR medial wall masks (Cell array of strings, MedialWallMasks{1} = path to LH, MedialWallMasks{2} = path to RH).
% HeadMesh: Path to tetrahedral head mesh (string).
% PositionUncertainity: Coil center positioning uncertainty (in mm, numeric). Optimal coil placement is defined as the position that maximizes the on target value on average within this distance. 
% CoilModel: Name of the coil model (string). Must be available in SimNIBS library ( /SimNIBS-3.2/simnibs/ccd-files/). 
% AbsoluteThreshold: This value (in V/m units) represents an assumed neural activation threshold (numeric). Default is 100 V/m. 
% DiDt: A range of stimulation intensities (the speed of variation of the current throughout the stimulating coil, in units of A/us, numeric). Default is DiDt = linspace(1,155,20) * 1e6, which corresponds *approximately* the possible range of realized 𝑑𝐼/𝑑𝑡 on our MagPro X100 machine when using the B70 coil.
% OutDir: Path to the output folder (string).
% Paths: Paths to folders that must be added to Matlab search path (cell array of strings). 

% define inputs 
Subject = 'ME01';
PialSurfs{1} = [DataDir '/' Subject '/anat/T1w/fsaverage_LR32k/' Subject '.L.pial.32k_fs_LR.surf.gii']; 
PialSurfs{2} = [DataDir '/' Subject '/anat/T1w/fsaverage_LR32k/' Subject '.R.pial.32k_fs_LR.surf.gii']; 
WhiteSurfs{1} = [DataDir '/' Subject '/anat/T1w/fsaverage_LR32k/' Subject '.L.white.32k_fs_LR.surf.gii'];
WhiteSurfs{2} = [DataDir '/' Subject '/anat/T1w/fsaverage_LR32k/' Subject '.R.white.32k_fs_LR.surf.gii']; 
MidthickSurfs{1} = [DataDir '/' Subject '/anat/T1w/fsaverage_LR32k/' Subject '.L.midthickness.32k_fs_LR.surf.gii']; 
MidthickSurfs{2} = [DataDir '/' Subject '/anat/T1w/fsaverage_LR32k/' Subject '.R.midthickness.32k_fs_LR.surf.gii']; 
MedialWallMasks{1} = [DataDir '/' Subject '/anat/MNINonLinear/fsaverage_LR32k/' Subject '.L.atlasroi.32k_fs_LR.shape.gii']; 
MedialWallMasks{2} = [DataDir '/' Subject '/anat/MNINonLinear/fsaverage_LR32k/' Subject '.R.atlasroi.32k_fs_LR.shape.gii']; 
VertexSurfaceArea = [DataDir '/' Subject '/anat/T1w/fsaverage_LR32k/' Subject '.midthickness_va.32k_fs_LR.dscalar.nii']; 
SkinSurf = [DataDir '/' Subject '/tans/HeadModel/m2m_' Subject '/Skin.surf.gii']; 
HeadMesh = [DataDir '/' Subject '/tans/HeadModel/' Subject '.msh']; 
OutDir = [DataDir '/' Subject '/tans/Network_Frontoparietal'];
PercentileThresholds = linspace(99.9,99,10);
CoilModel = 'MagVenture_MC_B70.nii.gz';  
PositionUncertainty = 5; 
DiDt = linspace(1,155,20) * 1e6; % A/us
AbsoluteThreshold = 100; % V/m

% isolate the target network again
TargetNetwork = FunctionalNetworks;
TargetNetwork.data(TargetNetwork.data~=9) = 0; % note: 9 == frontoparietal network map. 
TargetNetwork.data(TargetNetwork.data~=0) = 1; % binarize.

% run the "tans_optimize.m" module;
tans_optimize(Subject,TargetNetwork,[],PercentileThresholds,SkinVol,SkinSurf,SkinDistanceMatrix,VertexSurfaceArea,MidthickSurfs,WhiteSurfs,PialSurfs,MedialWallMasks,HeadMesh,PositionUncertainty,CoilModel,AbsoluteThreshold,DiDt,OutDir,Paths)

% generate an inverse transform (if needed); transform target coordinates to native T1w image space (which can be used for neuronav, for example)
system(['convert_xfm -omat ' DataDir '/' Subject '/anat/T1w/xfms/acpc2orig.mat -inverse ' DataDir '/' Subject '/anat/T1w/xfms/acpc.mat']);
system(['img2imgcoord -src ' DataDir '/' Subject '/anat/T1w/T1w_acpc.nii.gz -dest ' DataDir '/' Subject '/anat/T1w/T1w.nii.gz -xfm ' DataDir '/' Subject '/anat/T1w/xfms/acpc2orig.mat -mm ' OutDir '/Optimize/CoilCenterCoordinates.txt > ' OutDir '/Optimize/CoilCenterCordinatesNative.txt']);
system(['img2imgcoord -src ' DataDir '/' Subject '/anat/T1w/T1w_acpc.nii.gz -dest ' DataDir '/' Subject '/anat/T1w/T1w.nii.gz -xfm ' DataDir '/' Subject '/anat/T1w/xfms/acpc2orig.mat -mm ' OutDir '/Optimize/CoilOrientationCoordinates.txt > ' OutDir '/Optimize/CoilOrientationCordinatesNative.txt']);

%% visualize the extent to which the E-field "hotspot" is on target

% load the network labels and preset color information;
NetworkLabels = dataset('XLS',[DataDir '/' Subject '/pfm/' Subject '_FunctionalNetworks.xlsx']);
NetworkColors = [NetworkLabels.R NetworkLabels.G NetworkLabels.B]/255;
Networks = NetworkLabels.Network;

% generate the E-field hotspot associated with the optimal
E = ft_read_cifti_mod([DataDir '/' Subject '/tans/Network_Frontoparietal/Optimize/normE_BestCoilCenter+BestOrientation.dtseries.nii']);
[HotSpot] = tans_make_hotspot(E,PercentileThresholds);

% read in the surface vertex area information;
VertexSurfaceArea = ft_read_cifti_mod([DataDir '/' Subject '/anat/T1w/fsaverage_LR32k/' Subject '.midthickness_va.32k_fs_LR.dscalar.nii']);

% preallocate the on target variable;
OnTarget = zeros(length(PercentileThresholds),...
size(FunctionalNetworks.data,2));

% sweep the E-field thresholds
for t = 1:length(PercentileThresholds)
    
    % sweep all of 
    % the functional networks;
    for i = 1:length(Networks)
        ROI = FunctionalNetworks.data(1:59412)==i;
        OnTarget(t,i) = sum(VertexSurfaceArea.data(HotSpot(:,t+1) & ROI==1)) / sum(VertexSurfaceArea.data(HotSpot(:,t+1)==1));
    end
    
end

H = figure; % prellocate parent figure
set(H,'position',[1 1 400 250]); hold;

% plot the horizontal stacked bar graph;
h = barh([OnTarget ; nan(1,size(OnTarget,2)) ; mean(OnTarget,1)] * 100,1,'Stacked'); hold

% change bar colors
for ii = 1:length(Networks)
    h(ii).FaceColor = [NetworkColors(ii,1) NetworkColors(ii,2) NetworkColors(ii,3)];
end

% create some y-tick labels;
for i = 1:length(PercentileThresholds)
tmp{i} = [num2str(PercentileThresholds(i)) '%'];
end
tmp{length(PercentileThresholds)+1} = '';
tmp{length(PercentileThresholds)+2} = 'Average';

% make it "pretty";
yticks(1:length(PercentileThresholds)+2); xlim([0 100]);
ylim([0.5 length(PercentileThresholds)+2.5]);
xlabel('% of Total E-field Hotspot Surface Area');
set(gca,'FontName','Arial','FontSize',12,'TickLength',[0 0]);
yticklabels(tmp);

