function tans_simnibs(SearchGridCoords,HeadMesh,CoilModel,AngleResolution,SkinSurf,MidthickSurfs,WhiteSurfs,PialSurfs,MedialWallMasks,nThreads,OutDir,Paths)
% cjl; cjl2007@med.cornell.edu;

% The inputs to the function are:
% SearchGridCoords: Number of coil positions x 3 numeric array generated by tans_searchgrid function. 
% HeadMesh: Path to tetrahedral head mesh (string).
% CoilModel: Name of the coil model (string). Must be available in SimNIBS library ( /SimNIBS-3.2/simnibs/ccd-files/). 
% AngleResolution: Inter angle resolution used in the search grid (numeric).
% SkinSurf: Path to the skin surface geometry file (string).
% MidthickSurfs: Paths to low (32k) dimensional FS_LR midthickness surfaces (Cell array of strings, MidthickSurfs{1} = path to LH, MidthickSurfs{2} = path to RH).
% WhiteSurfs: Paths to low (32k) dimensional FS_LR white surfaces (Cell array of strings, WhiteSurfs{1} = path to LH, WhiteSurfs{2} = path to RH).
% PialSurfs: Paths to low (32k) dimensional FS_LR pial surfaces (Cell array of strings, PialSurfs{1} = path to LH, PialSurfs{2} = path to RH).
% MedialWallMasks: Paths to low (32k) dimensional FS_LR medial wall masks (Cell array of strings, MedialWallMasks{1} = path to LH, MedialWallMasks{2} = path to RH).
% nThreads: The number of parallel pools that will be used in Matlab (numeric).
% OutDir: Path to the output folder (string).
% Paths: Paths to folders that must be added to Matlab search path (cell array of strings). 

% define some directories;
addpath(genpath(Paths{1})); % define the path to SimNibs software
addpath(genpath(Paths{2})); % define the path to the folder containing "ft_read_cifti" / "gifti" functions

rng(44); % for reproducibility;
warning ('off','all'); % turn off annoying warnings;

% load skin mesh generated
% by "tans_headmodels.m"
SkinSurf = gifti(SkinSurf); 

% define a parpool;
pool = parpool('local',nThreads);

% sweep the search space;
parfor i = 1:size(SearchGridCoords,1)
    
    try % note: sometimes the scalp position is weird and the simulation will fail. 
        % So we have a try/catch here to prevent the whole thing from failing.
    
    % Initialize a session
    s = sim_struct('SESSION');
    
    % Name of head mesh
    s.fnamehead = HeadMesh;
    
    % Output folder
    s.pathfem = [OutDir '/SearchGrid/Simulation_' sprintf('%05d',i) '/'];
    
    % Initialize a list of TMS simulations
    s.poslist{1} = sim_struct('TMSLIST');
    
    % specify the coil model;
    s.poslist{1}.fnamecoil = CoilModel;
    
    % generate an approximate circle around the center position (this is a bit crude but works fine); 
    D = pdist2(SkinSurf.vertices,[SearchGridCoords(i,1), SearchGridCoords(i,2), SearchGridCoords(i,3)]);
    A = find(D < 19); % inner diameter;
    B = find(D < 20); % outer diameter;
    Circle = SkinSurf.vertices(B(~ismember(B,A)),:); % find the difference;
    
    % cluster 3d coordinates;
    [~,yDirs] = kmeans(Circle,round(360 / AngleResolution));
    
    % sort coordinates in a radial fashion;
    [yDirs] = SortCircleCoords(yDirs); % 
    yDirs = yDirs(1:ceil(size(yDirs,1)/2),:); % sample one half of the circle; this reduces the number of simulations you need to run, while still capturing the max on-target value for this coil center. 

    % sweep a range of 
    % coil orientations;
    for ii = 1:length(yDirs)
    
    % specify & save coil position;
    s.poslist{1}.pos(ii).centre = [SearchGridCoords(i,1), SearchGridCoords(i,2), SearchGridCoords(i,3)];
    s.poslist{1}.pos(ii).pos_ydir = [yDirs(ii,1), yDirs(ii,2), yDirs(ii,3)];
    s.poslist{1}.pos(ii).didt = 1 * 1e6; % A/us

    end
    
    % write to volume;
    s.map_to_vol = true;
    s.fields = 'e'; % normE only;
    
    % run the
    % simulation;
    run_simnibs(s);
    
    % merge all the volumes into a single 4D file;
    system(['fslmerge -t ' OutDir '/SearchGrid/Simulation_' sprintf('%05d',i) '/subject_volumes/normE.nii.gz ' OutDir '/SearchGrid/Simulation_' sprintf('%05d',i) '/subject_volumes/*normE.nii.gz']);
    
    % map concatenated volume to the 32k surface;
    system(['wb_command -volume-to-surface-mapping ' OutDir '/SearchGrid/Simulation_' sprintf('%05d',i) '/subject_volumes/normE.nii.gz  ' MidthickSurfs{1} ' ' OutDir '/SearchGrid/Simulation_' sprintf('%05d',i) '/subject_volumes/normE.L.32k_fs_LR.shape.gii -ribbon-constrained ' WhiteSurfs{1} ' ' PialSurfs{1}]); % note: on 1/11/22 CJL removed -interpolate ENCLOSING at the end of this command. 
    system(['wb_command -volume-to-surface-mapping ' OutDir '/SearchGrid/Simulation_' sprintf('%05d',i) '/subject_volumes/normE.nii.gz  ' MidthickSurfs{2} ' ' OutDir '/SearchGrid/Simulation_' sprintf('%05d',i) '/subject_volumes/normE.R.32k_fs_LR.shape.gii -ribbon-constrained ' WhiteSurfs{2} ' ' PialSurfs{2}]);
    system(['wb_command -metric-mask ' OutDir '/SearchGrid/Simulation_' sprintf('%05d',i) '/subject_volumes/normE.L.32k_fs_LR.shape.gii ' MedialWallMasks{1} ' ' OutDir '/SearchGrid/Simulation_' sprintf('%05d',i) '/subject_volumes/normE.L.32k_fs_LR.shape.gii']); 
    system(['wb_command -metric-mask ' OutDir '/SearchGrid/Simulation_' sprintf('%05d',i) '/subject_volumes/normE.R.32k_fs_LR.shape.gii ' MedialWallMasks{2} ' ' OutDir '/SearchGrid/Simulation_' sprintf('%05d',i) '/subject_volumes/normE.R.32k_fs_LR.shape.gii']);
    system(['wb_command -cifti-create-dense-timeseries ' OutDir '/SearchGrid/Simulation_' sprintf('%05d',i) '/subject_volumes/normE.dtseries.nii -left-metric ' OutDir '/SearchGrid/Simulation_' sprintf('%05d',i) '/subject_volumes/normE.L.32k_fs_LR.shape.gii -roi-left ' MedialWallMasks{1} ' -right-metric ' OutDir '/SearchGrid/Simulation_' sprintf('%05d',i) '/subject_volumes/normE.R.32k_fs_LR.shape.gii -roi-right ' MedialWallMasks{2}]);
    system(['mv ' OutDir '/SearchGrid/Simulation_' sprintf('%05d',i) '/subject_volumes/normE.dtseries.nii ' OutDir '/SearchGrid/normE_' sprintf('%05d',i) '.dtseries.nii']);
    system(['mv ' OutDir '/SearchGrid/Simulation_' sprintf('%05d',i) '/*.log ' OutDir '/SearchGrid/Simulation_' sprintf('%05d',i) '.log']);

    % remove some large intermediate files to save disk space;
    system(['rm -rf ' OutDir '/SearchGrid/Simulation_' sprintf('%05d',i) '/']);
  
    catch
    end
    
end

% delete the
% parpool;
delete(pool);

end

% sub-functions
function [Output] = SortCircleCoords(Coords)
% cjl; purpose: take a set of 3d coordinates 
% and order them in a radial fashion;

% preallocate;
Output = zeros(size(Coords));

% preallocate;
Used = [];
Idx = 1;

% sweep the coordinates;
for i = 1:size(Coords,1)
    D = pdist2(Coords(Idx,:),Coords);
    D(Used) = nan;
    Idx = find(D==min(nonzeros(D)));
    Idx = Idx(1); % in case there is more than one;
    Output(i,:) = Coords(Idx,:);
    Used = [Used Idx];
end

end



