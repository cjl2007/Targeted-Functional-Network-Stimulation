%% Example use of Targeted Functional Network Stimulation ("TANS")
% "Automated optimization of TMS coil placement for personalized functional
% network engagement" - Lynch et al., 2022 (Neuron)

% define some paths
Paths{1} = '/SimNIBS-4.0'; % download from https://simnibs.github.io/simnibs/build/html/index.html
Paths{2} = '/Utilities/msc'; % this folder contains ft_read / gifti functions for reading and writing cifti files (e.g., https://github.com/MidnightScanClub/MSCcodebase).
Paths{3} = '/Targeted-Functional-Network-Stimulation-main/'; %

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
[status,~] = system('charm -version'); % simnibs: note charm program replaced headreco in SimNIBS v4.0.
[status,~] = system('wb_command -version'); % connectome workbench
[status,~] = system('flirt -version'); % fsl

% If successful, each of the commands below should return status as 2.

% Confirm that functions for reading and writing
% CIFTI and GIFTI files are also available
status = exist('ft_read_cifti_mod','file');
status = exist('gifti','file');

% define the path to directory containing example data;
DataDir = '/home/charleslynch/Targeted-Functional-Network-Stimulation-main/Data/'; % note: should be obvious, but user needs to manually define this

% sweep
% subjects
for i = 1:6
    
    % Create the head model
    % Timing: Typically 1-2 hours.
    
    % Use the tans_headmodels function to create a head model for electric field (E-field) modeling
    % using the CHARM method (Puonti et al. 2020). In addition, this function generates skin
    % surface geometry files (*.surf.gii) for visualizing certain results.
    
    % The inputs to the function are:
    % Subject: The subject ID (string).
    % T1w: Path to T1w-weighted anatomical image (string).
    % T2w: Path to T2w-weighted anatomical image (string).
    % OutDir : Path to the output folder (string).
    % Paths: Paths to folders that must be added to Matlab search path (cell array of strings).
    
    % define inputs
    Subject = ['ME0' num2str(i)];
    T1w = [DataDir '/' Subject '/anat/T1w/T1w_acpc_dc_restore.nii.gz'];
    T2w = [DataDir '/' Subject '/anat/T1w/T2w_acpc_dc_restore.nii.gz'];
    OutDir = [DataDir '/' Subject '/tans'];
    
    % run the tans_headmodels function;
    tans_headmodels(Subject,T1w,T2w,OutDir,Paths);
    
    % Identify the target network patch.
    % Timing: Typically less than 1 minute.
    
    % Functional networks usually consist of many nodes spanning multiple cortical zones in both hemispheres.
    % It is unclear a priori which part of the network TMS should be delivered to maximize stimulation
    % specificity. Use the tans_roi function to find the largest patch of the target network on a
    % gyral crown within a search space specified by the user (for example, lateral prefrontal cortex).
    
    % The inputs to the function are:
    % TargetNetwork: A CIFTI file containing the functional network of interest (structure array). Non-zero values in TargetNetwork.data are considered target network vertices.
    % MidthickSurfs: Paths to low (32k) dimensional FS_LR midthickness surfaces (Cell array of strings, MidthickSurfs{1} = path to LH, MidthickSurfs{2} = path to RH).
    % VertexSurfaceArea: A CIFTI file containing the vertex surface areas (structure array).
    % Sulc:  CIFTI file containing sulcal depth information (structure array).
    % SearchSpace: A CIFTI file containing a binarized mask representing the search space (structure array).
    % OutDir: Path to the output folder (string).
    % Paths: Paths to folders that must be added to Matlab search path (cell array of strings).
    
    % read in the functional networks for this individual;
    FunctionalNetworks = ft_read_cifti_mod([DataDir '/' Subject '/pfm/'...
    Subject '_FunctionalNetworks.dtseries.nii']);
    
    % isolate the target network
    TargetNetwork = FunctionalNetworks;
    TargetNetwork.data(TargetNetwork.data~=9) = 0; % note: 9 == frontoparietal network map.
    TargetNetwork.data(TargetNetwork.data~=0) = 1; % binarize.
    
    % load cortical mask representing the lateral PFC;
    SearchSpace = ft_read_cifti_mod('LPFC_LH+RH.dtseries.nii');
    
    % load sulcal depth information;
    Sulc = ft_read_cifti_mod([DataDir '/' Subject '/anat/MNINonLinear/fsaverage_LR32k/' Subject '.sulc.32k_fs_LR.dscalar.nii']);
    BrainStructure = SearchSpace.brainstructure; % extract the brain structure index
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
    
    % Make a search grid on the scalp above the target network patch centroid.
    % Timing: Typically 5-10 minutes.
    
    % Use the tans_searchgrid function to generate an array of three-dimensional coordinates
    % representing a search grid on the scalp directly above the centroid of the target network patch.
    
    % The inputs to the function are:
    % TargetPatch: A CIFTI file containing the target network patch (structure array). Non-zero values in TargetPatch.data are considered target network patch vertices.
    % PialSurfs: Paths to low (32k) dimensional FS_LR pial surfaces (Cell array of strings, PialSurfs{1} = path to LH, PialSurfs{2} = path to RH).
    % SkinSurf: Path to the skin surface geometry file (string).
    % SearchGridRadius: Radius (in mm) of the search grid on the scalp directly above the centroid of the target network patch.
    % GridSpacing: The (approximate) distance between pairs of vertices in the search grid, in mm.
    % OutDir: Path to the output folder (string).
    % Paths: Paths to folders that must be added to Matlab search path (cell array of strings).
    
    % define inputs
    OutDir = [DataDir '/' Subject '/tans/Network_Frontoparietal'];
    TargetNetworkPatch = ft_read_cifti_mod([OutDir '/ROI/TargetNetworkPatch.dtseries.nii']);
    PialSurfs{1} = [DataDir '/' Subject '/anat/T1w/fsaverage_LR32k/' Subject '.L.pial.32k_fs_LR.surf.gii'];
    PialSurfs{2} = [DataDir '/' Subject '/anat/T1w/fsaverage_LR32k/' Subject '.R.pial.32k_fs_LR.surf.gii'];
    SkinSurf = [DataDir '/' Subject '/tans/HeadModel/m2m_' Subject '/Skin.surf.gii'];
    SearchGridRadius = 20;
    GridSpacing = 2;
    
    % run the tans_searchgrid function
    [SubSampledSearchGrid,FullSearchGrid] = tans_searchgrid(TargetNetworkPatch,PialSurfs,SkinSurf,GridSpacing,SearchGridRadius,OutDir,Paths);
    
    % Perform electric field modeling iteratively in the search grid.
    % Timing: Depends on the total number of simulations performed and available parallel pools in Matlab, typically several hours.
    
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
    HeadMesh = [DataDir '/' Subject '/tans/HeadModel/m2m_' Subject '/' Subject '.msh'];
    CoilModel = [Paths{1} '/simnibs_env/lib/python3.9/site-packages/simnibs/resources/coil_models/Drakaki_BrainStim_2022/MagVenture_MCF-B65.ccd']; % note: SimNIBS 4.0 seems to require full path, whereas earlier versions did not.
    DistanceToScalp = 2;
    AngleResolution = 30;
    nThreads = 20;
    
    % run the tans_simnibs function
    tans_simnibs(SearchGridCoords,HeadMesh,CoilModel,AngleResolution,DistanceToScalp,SkinSurf,MidthickSurfs,WhiteSurfs,PialSurfs,MedialWallMasks,nThreads,OutDir,Paths);
    
    % Note: A fixed stimulation intensity (𝑑𝐼/𝑑𝑡 = 1 A/µs) is used during this stage of TANS because the strength of the E-field
    % varies linearly with 𝑑𝐼/𝑑𝑡 (the speed of variation of the current throughout the coil) and, for this reason,
    % has no effect on its spatial distribution (including where it is maximal relative to the target).
    
    % Find the coil placement that best aligns the electric field hotspot with the target network
    % Timing: Depends on the total number of simulations performed, typically 1-2 hours.
    
    % Use the tans_optimize function to find which of the coil placements submitted to the tans_simnibs function
    % produced the greatest on-target value (calculated as the proportion of the E-field hotspot inside the
    % target network). In other words, which coil placement was best for maximizing stimulation specificity.
    
    % The inputs to the function are:
    % Subject: The subject ID (string).
    % TargetNetwork: A CIFTI file containing the functional network of interest (structure array). Non-zero values in TargetNetwork.data are considered target network vertices.
    % AvoidanceRegion: A CIFTI file indexing the brain regions or networks of no interest (structure array). Non-zero values in AvoidanceRegion.data are considered non-target network vertices.
    % PercentileThresholds: A range of percentiles used for operationalizing the E-field hotspot (numeric). For example, linspace(99.9,99,10).
    % SkinSurf: Path to the skin surface geometry file (string).
    % VertexSurfaceArea: A CIFTI file containing the vertex surface areas (structure array).
    % MidthickSurfs: Paths to low (32k) dimensional FS_LR midthickness surfaces (Cell array of strings, MidthickSurfs{1} = path to LH, MidthickSurfs{2} = path to RH).
    % WhiteSurfs: Paths to low (32k) dimensional FS_LR white surfaces (Cell array of strings, WhiteSurfs{1} = path to LH, WhiteSurfs{2} = path to RH).
    % PialSurfs: Paths to low (32k) dimensional FS_LR pial surfaces (Cell array of strings, PialSurfs{1} = path to LH, PialSurfs{2} = path to RH).
    % MedialWallMasks: Paths to low (32k) dimensional FS_LR medial wall masks (Cell array of strings, MedialWallMasks{1} = path to LH, MedialWallMasks{2} = path to RH).
    % AngleResolution: Inter angle resolution used to fine tune coil orientation (numeric).
    % HeadMesh: Path to tetrahedral head mesh (string).
    % PositionUncertainity: Coil center positioning uncertainty (in mm, numeric). Optimal coil placement is defined as the position that maximizes the on target value on average within this distance.
    % CoilModel: Path to the coil model (string). Must be available in SimNIBS library (/SimNIBS-3.2/simnibs/ccd-files/).
    % OutDir: Path to the output folder (string).
    % Paths: Paths to folders that must be added to Matlab search path (cell array of strings).
    
    % define inputs
    PialSurfs{1} = [DataDir '/' Subject '/anat/T1w/fsaverage_LR32k/' Subject '.L.pial.32k_fs_LR.surf.gii'];
    PialSurfs{2} = [DataDir '/' Subject '/anat/T1w/fsaverage_LR32k/' Subject '.R.pial.32k_fs_LR.surf.gii'];
    WhiteSurfs{1} = [DataDir '/' Subject '/anat/T1w/fsaverage_LR32k/' Subject '.L.white.32k_fs_LR.surf.gii'];
    WhiteSurfs{2} = [DataDir '/' Subject '/anat/T1w/fsaverage_LR32k/' Subject '.R.white.32k_fs_LR.surf.gii'];
    MidthickSurfs{1} = [DataDir '/' Subject '/anat/T1w/fsaverage_LR32k/' Subject '.L.midthickness.32k_fs_LR.surf.gii'];
    MidthickSurfs{2} = [DataDir '/' Subject '/anat/T1w/fsaverage_LR32k/' Subject '.R.midthickness.32k_fs_LR.surf.gii'];
    VertexSurfaceArea = ft_read_cifti_mod([DataDir '/' Subject '/anat/T1w/fsaverage_LR32k/' Subject '.midthickness_va.32k_fs_LR.dscalar.nii']);
    MedialWallMasks{1} = [DataDir '/' Subject '/anat/MNINonLinear/fsaverage_LR32k/' Subject '.L.atlasroi.32k_fs_LR.shape.gii'];
    MedialWallMasks{2} = [DataDir '/' Subject '/anat/MNINonLinear/fsaverage_LR32k/' Subject '.R.atlasroi.32k_fs_LR.shape.gii'];
    SearchGrid = [DataDir '/' Subject '/tans/Network_Frontoparietal/SearchGrid/SubSampledSearchGrid.shape.gii'];
    SkinFile = [DataDir '/' Subject '/tans/HeadModel/m2m_' Subject '/Skin.surf.gii'];
    HeadMesh = [DataDir '/' Subject '/tans/HeadModel/m2m_' Subject '/' Subject '.msh'];
    OutDir = [DataDir '/' Subject '/tans/Network_Frontoparietal/'];
    PercentileThresholds = linspace(99.9,99,10);
    CoilModel = [Paths{1} '/simnibs_env/lib/python3.9/site-packages/simnibs/resources/coil_models/Drakaki_BrainStim_2022/MagVenture_MCF-B65.ccd']; % note: SimNIBS 4.0 seems to require full path, whereas earlier versions did not.
    DistanceToScalp = 2;
    Uncertainty = 5;
    AngleResolution = 5;
    
    % isolate the target network again
    TargetNetwork = FunctionalNetworks;
    TargetNetwork.data(TargetNetwork.data~=9) = 0; % note: 9 == frontoparietal network map.
    TargetNetwork.data(TargetNetwork.data~=0) = 1; % binarize.
    
    % run the "tans_optimize.m" module;
    tans_optimize(Subject,TargetNetwork,[],PercentileThresholds,SearchGrid,DistanceToScalp,SkinSurf,VertexSurfaceArea,MidthickSurfs,WhiteSurfs,PialSurfs,MedialWallMasks,HeadMesh,AngleResolution,Uncertainty,CoilModel,OutDir,Paths);
    
    % Estimate the stimulation intensity level that maximizes stimulation
    % specificity after some assumptions regarding neural activation threshold and
    % a min. e-field hotspot size.
    % Timing: 1 minute.
    
    % define some new input variables;
    NetworkLabels = [DataDir '/' Subject '/pfm/' Subject '_FunctionalNetworks.xlsx'];
    magnE = [OutDir '/Optimize/magnE_BestCoilCenter+BestOrientation.dtseries.nii'];
    DiDt = linspace(1,155,20) * 1e6; % A/us
    AbsoluteThreshold = 100; % V/m
    MinHotSpotSize = 1000; % surface area (mm2)
    
    % run the "tans_dose.m" module
    [OnTarget,Penalty,HotSpotSize] =tans_dose(magnE,VertexSurfaceArea,DiDt,AbsoluteThreshold,MinHotSpotSize,TargetNetwork,[],FunctionalNetworks,NetworkLabels,[OutDir '/Optimize/'],Paths);
    
end
