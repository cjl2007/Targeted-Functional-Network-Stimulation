function tans_optimize(Subject,TargetNetwork,AvoidanceRegion,PercentileThresholds,SearchGrid,DistanceToScalp,SkinFile,VertexSurfaceArea,MidthickSurfs,WhiteSurfs,PialSurfs,MedialWallMasks,HeadMesh,AngleResolution,Uncertainty,CoilModel,OutDir,Paths)
% cjl; cjl2007@med.cornell.edu;
%
% The inputs to the function are:
% Subject: The subject ID (string).
% TargetNetwork: A CIFTI file containing the functional network of interest (structure array). Non-zero values in TargetNetwork.data are considered target network vertices.
% AvoidanceRegion: A CIFTI file indexing the brain regions or networks of no interest (structure array). Non-zero values in AvoidanceRegion.data are considered non-target network vertices. 
% PercentileThresholds: A range of percentiles used for operationalizing the E-field hotspot (numeric). For example, linspace(99.9,99,10).
% SearchGrid: Path to search grid gifti file. 
% DistanceToScalp: Distance between coil center and scalp surface.
% SkinFile: Path to the skin surface geometry file (string). 
% VertexSurfaceArea: A CIFTI file containing the vertex surface areas (structure array). 
% MidthickSurfs: Paths to low (32k) dimensional FS_LR midthickness surfaces (Cell array of strings, MidthickSurfs{1} = path to LH, MidthickSurfs{2} = path to RH).
% WhiteSurfs: Paths to low (32k) dimensional FS_LR white surfaces (Cell array of strings, WhiteSurfs{1} = path to LH, WhiteSurfs{2} = path to RH).
% PialSurfs: Paths to low (32k) dimensional FS_LR pial surfaces (Cell array of strings, PialSurfs{1} = path to LH, PialSurfs{2} = path to RH).
% MedialWallMasks: Paths to low (32k) dimensional FS_LR medial wall masks (Cell array of strings, MedialWallMasks{1} = path to LH, MedialWallMasks{2} = path to RH).
% AngleResolution: Inter angle resolution for fine tuning coil orientation (numeric).
% HeadMesh: Path to tetrahedral head mesh (string).
% Uncertainty: Coil center positioning uncertainty (in mm, numeric). Optimal coil placement is defined as the position that maximizes the on target value on average within this distance. 
% CoilModel: Path to the coil model (string). 
% OutDir: Path to the output folder (string).
% Paths: Paths to folders that must be added to Matlab search path (cell array of strings). 

% add some directories 
% to the search path;
for i = 1:length(Paths)
addpath(genpath(Paths{i})); % 
end

rng(44); % for reproducibility;
warning ('off','all'); % turn off annoying warnings; users could comment this line out if they want. 
SmoothingFactor = 0.85; % might need to increase if the search grid is very sparse.

% make the optimize dir.
mkdir([OutDir '/Optimize']);

% SkinFile == path to skin surface file;
% Skin == loaded version
Skin = gifti(SkinFile); 

% evaluate E-fields associated with all the coil placements;

% load the search 
% grid metric file;
G = gifti(SearchGrid); 
SearchGridVertices = find(G.cdata~=0); % 

% count the number of files (note: some attempts will have failed,
% for whatever reason, so we cant assume that file X is simulation X;
files = dir([OutDir '/SearchGrid/magnE_*.dtseries.nii']);

% read the first file & count how many orientations were attempted
MagnE = ft_read_cifti_mod([OutDir '/SearchGrid/' files(1).name]);
nCols = size(MagnE.data,2);

% preallocate;
OnTarget = zeros(length(SearchGridVertices),nCols,length(PercentileThresholds)); % "OnTaget" variable (% of E-field hotspot that contains target network vertices);
Penalty = zeros(length(SearchGridVertices),nCols,length(PercentileThresholds)); % "Penalty" variable (% of E-field hotspot that contains avoidance region / network vertices);

% if no avoidance 
% region is specified; 
if isempty(AvoidanceRegion)
AvoidanceRegion = TargetNetwork; % preallocate
AvoidanceRegion.data = zeros(size(AvoidanceRegion.data,1));
end

% sweep the search space;
for i = 1:length(files)
   
    % read in the CIFTI file;
    MagnE = ft_read_cifti_mod([OutDir '/SearchGrid/' files(i).name]);
    system(['rm ' OutDir '/SearchGrid/' files(i).name]); % remove intermediate files;

    % make sure we have the correct point in the space grid;
    tmp = strsplit(files(i).name,{'_','.'}); % note: sometimes a given point in the search grid will fail;
    idx = str2double(tmp{2});
    
    % sweep the coil orientations;
    for ii = 1:size(MagnE.data,2)
        
        % sweep all of the thresholds;
        for iii = 1:length(PercentileThresholds)
            
            HotSpot = MagnE.data(1:59412,ii) > prctile(MagnE.data(1:59412,ii),PercentileThresholds(iii)); % this is the hotspot
            OnTarget(idx,ii,iii) = (sum(VertexSurfaceArea.data(HotSpot&TargetNetwork.data(1:59412,1)==1)) / sum(VertexSurfaceArea.data(HotSpot))); 
            Penalty(idx,ii,iii) = (sum(VertexSurfaceArea.data(HotSpot&AvoidanceRegion.data(1:59412,1)==1)) / sum(VertexSurfaceArea.data(HotSpot)));
            
        end
        
    end

    % clear 
    % variable
    clear MagnE
    
end

% save some variables;
save([OutDir '/Optimize/CoilCenter_OnTarget'],'OnTarget');
save([OutDir '/Optimize/CoilCenter_Penalty'],'Penalty');

% average accross the e-field thresholds;
AvgOnTarget = mean(OnTarget,3); % on-target;
AvgPenalty = mean(Penalty,3); % penalty;
AvgPenalizedOnTarget = mean(OnTarget,3) - mean(Penalty,3); % on-target - penalty; relative on-target value used for optimization

G_OnTarget = G; % preallocate;
G_OnTarget.cdata = zeros(size(G_OnTarget.cdata)); % blank slate;

% write out the on-target metric file;
G_OnTarget.cdata(SearchGridVertices) = max(AvgOnTarget,[],2); % average across orientations, for now.
save(G_OnTarget,[OutDir '/Optimize/CoilCenter_OnTarget.shape.gii']); % write out the on-target metric file;
system(['wb_command -metric-smoothing ' SkinFile ' ' OutDir '/Optimize/CoilCenter_OnTarget.shape.gii ' num2str(SmoothingFactor) ' ' OutDir '/Optimize/CoilCenter_OnTarget_s' num2str(SmoothingFactor) '.shape.gii -fix-zeros']);
G_OnTarget = gifti([OutDir '/Optimize/CoilCenter_OnTarget_s' num2str(SmoothingFactor) '.shape.gii']); % read in the smoothed file;

% write out the penalty metric file;
G_Penalty = G_OnTarget; % preallocate
G_Penalty.cdata = zeros(size(G_OnTarget.cdata)); % blank slate;
G_Penalty.cdata(SearchGridVertices) = max(AvgPenalty,[],2); % average across orientations, for now.
save(G_Penalty,[OutDir '/Optimize/CoilCenter_Penalty.shape.gii']); % write out the on-target metric file;
system(['wb_command -metric-smoothing ' SkinFile ' ' OutDir '/Optimize/CoilCenter_Penalty.shape.gii ' num2str(SmoothingFactor) ' ' OutDir '/Optimize/CoilCenter_Penalty_s' num2str(SmoothingFactor) '.shape.gii -fix-zeros']);

% write out the penalized on-target metric file;
G_PenalizedOnTarget = G_OnTarget; % preallocate
G_PenalizedOnTarget.cdata = zeros(size(G_OnTarget.cdata)); % blank slate;
G_PenalizedOnTarget.cdata(SearchGridVertices) = max(AvgPenalizedOnTarget,[],2); % average across orientations, for now.
save(G_PenalizedOnTarget,[OutDir '/Optimize/CoilCenter_PenalizedOnTarget.shape.gii']); % write out the relative on-target metric file;
system(['wb_command -metric-smoothing ' SkinFile ' ' OutDir '/Optimize/CoilCenter_PenalizedOnTarget.shape.gii ' num2str(SmoothingFactor) ' ' OutDir '/Optimize/CoilCenter_PenalizedOnTarget_s' num2str(SmoothingFactor) '.shape.gii -fix-zeros']);
G_PenalizedOnTarget = gifti([OutDir '/Optimize/CoilCenter_PenalizedOnTarget_s' num2str(SmoothingFactor) '.shape.gii']); % read in the smoothed file;

% preallocate the overall quality (adjusted for some amoutn of error
% anticipated from neuronavigation imprecision).
PenalizedOnTarget_ErrorAdjusted = zeros(length(SearchGridVertices),1);

% sweep the search grid vertices;
for i = 1:length(SearchGridVertices)
    D = pdist2(Skin.vertices,Skin.vertices(SearchGridVertices(i),:));
    PenalizedOnTarget_ErrorAdjusted(i) = mean(G_PenalizedOnTarget.cdata(D <= Uncertainty)); % average value of all vertices within the specified distance of vertex "i" 
end

% define the best coil center placement;
Idx = find(PenalizedOnTarget_ErrorAdjusted==max(PenalizedOnTarget_ErrorAdjusted)); 
CoilCenterVertex = SearchGridVertices(Idx(1)); % this is the final coil placement site;
CoilCenterCoords = Skin.vertices(CoilCenterVertex,:);
system(['echo ' num2str(CoilCenterCoords(1)) ' ' num2str(CoilCenterCoords(2)) ' ' num2str(CoilCenterCoords(3)) ' > ' OutDir '/Optimize/CoilCenterCoordinates.txt']); % write coordinates out to .txt file;

% make a foci on the skin surface;
system(['echo Target > ' OutDir '/Optimize/tmp.txt']);
system(['echo 1 0 0 ' num2str(round(CoilCenterCoords(1))) ' ' num2str(round(CoilCenterCoords(2))) ' ' num2str(round(CoilCenterCoords(3))) ' >> ' OutDir '/Optimize/tmp.txt']);
system(['wb_command -foci-create  ' OutDir '/Optimize/CoilCenter.foci -class CoilCenter ' OutDir '/Optimize/tmp.txt ' SkinFile]);
system(['rm ' OutDir '/Optimize/tmp.txt']); % remove intermediate file;

%% fine-tune the coil orientation

% Initialize a session
s = sim_struct('SESSION');

% Name of head mesh
s.fnamehead = HeadMesh;

% Output folder
s.pathfem = [OutDir '/Optimize/Simulation/'];

% Initialize a list of TMS simulations
s.poslist{1} = sim_struct('TMSLIST');

% Select coil
s.poslist{1}.fnamecoil = CoilModel;

% specify coil positioning
s.poslist{1}.pos(1).centre = CoilCenterCoords;

% generate an approximate circle around the center position;
D = pdist2(Skin.vertices,CoilCenterCoords);
A = find(D < 19); % inner diameter;
B = find(D < 20); % outer diameter;
Circle = Skin.vertices(B(~ismember(B,A)),:); % find the difference;

% convert to a single ring (remove any "pile up");
[~,Circle] = kmeans(Circle,round(360 / AngleResolution));

% sort coordinates in a radial fashion & average across neighboring 
% positions to achieve more regular spacing between reference points;
[RefDirs] = SortCircleCoords(Circle,ceil(size(Circle,1) * 0.15)); %

% preallocate;
RefDirsVertices = ...
zeros(size(RefDirs,1),1);

% sweep a range of
% coil orientations;
for i = 1:length(RefDirs)
    
    % specify & save coil position;
    s.poslist{1}.pos(i).centre = CoilCenterCoords;
    s.poslist{1}.pos(i).pos_ydir = [RefDirs(i,1), RefDirs(i,2), RefDirs(i,3)];
    s.poslist{1}.pos(i).distance = DistanceToScalp; % 

    % while we are at it, lets log 
    % the corresponding skin vertex numbers;
    Tmp = pdist2(Skin.vertices,RefDirs(i,:));
    Idx = find(Tmp==min(Tmp));
    RefDirsVertices(i) = Idx(1);
    
end
    
% write to volume;
s.map_to_vol = true;
s.fields = 'e'; % magnE only;

% remove the dir. if it exists;
if exist([OutDir '/Optimize/Simulation/'],'dir')
system(['rm -rf ' OutDir '/Optimize/Simulation/']);
end

%run simulation;
run_simnibs(s);

% map concatenated volume to the 32k surface;
system(['fslmerge -t ' OutDir '/Optimize/Simulation/subject_volumes/magnE.nii.gz ' OutDir '/Optimize/Simulation/subject_volumes/' Subject '*_magnE.nii.gz']);
system(['rm ' OutDir '/Optimize/Simulation/subject_volumes/' Subject '*_magnE.nii.gz']); % remove intermediate files;
system(['wb_command -volume-to-surface-mapping ' OutDir '/Optimize/Simulation/subject_volumes/magnE.nii.gz  ' MidthickSurfs{1} ' ' OutDir '/Optimize/Simulation/subject_volumes/magnE.L.32k_fs_LR.shape.gii -ribbon-constrained ' WhiteSurfs{1} ' ' PialSurfs{1} ' -interpolate ENCLOSING_VOXEL']);
system(['wb_command -volume-to-surface-mapping ' OutDir '/Optimize/Simulation/subject_volumes/magnE.nii.gz  ' MidthickSurfs{2} ' ' OutDir '/Optimize/Simulation/subject_volumes/magnE.R.32k_fs_LR.shape.gii -ribbon-constrained ' WhiteSurfs{2} ' ' PialSurfs{2} ' -interpolate ENCLOSING_VOXEL']);
system(['wb_command -metric-mask ' OutDir '/Optimize/Simulation/subject_volumes/magnE.L.32k_fs_LR.shape.gii ' MedialWallMasks{1} ' ' OutDir '/Optimize/Simulation/subject_volumes/magnE.L.32k_fs_LR.shape.gii']);
system(['wb_command -metric-mask ' OutDir '/Optimize/Simulation/subject_volumes/magnE.R.32k_fs_LR.shape.gii ' MedialWallMasks{2} ' ' OutDir '/Optimize/Simulation/subject_volumes/magnE.R.32k_fs_LR.shape.gii']);
system(['wb_command -cifti-create-dense-timeseries ' OutDir '/Optimize/Simulation/subject_volumes/magnE.dtseries.nii -left-metric ' OutDir '/Optimize/Simulation/subject_volumes/magnE.L.32k_fs_LR.shape.gii -roi-left ' MedialWallMasks{1} ' -right-metric ' OutDir '/Optimize/Simulation/subject_volumes/magnE.R.32k_fs_LR.shape.gii -roi-right ' MedialWallMasks{2}]);
system(['rm ' OutDir '/Optimize/Simulation/subject_volumes/*shape*']);

% rename the cifti file and remove some intermediate files;
system(['mv ' OutDir '/Optimize/Simulation/subject_volumes/magnE.dtseries.nii '...
OutDir '/Optimize/magnE_BestCoilCenter+AllOrientations.dtseries.nii']);

% read the cifti file & count how many orientations were attempted
MagnE = ft_read_cifti_mod([OutDir '/Optimize/magnE_BestCoilCenter+AllOrientations.dtseries.nii']);
nCols = size(MagnE.data,2);

% preallocate;
OnTarget = zeros(1,nCols,length(PercentileThresholds)); % "OnTaget" variable (% of E-field hotspot that contains target network vertices);
Penalty = zeros(1,nCols,length(PercentileThresholds)); % "Penalty" variable (% of E-field hotspot that contains avoidance region / network vertices);

% sweep the coil orientations;
for i = 1:size(MagnE.data,2)
    
    % sweep all of the thresholds;
    for ii = 1:length(PercentileThresholds)

        HotSpot = MagnE.data(1:59412,i) > prctile(MagnE.data(1:59412,i),PercentileThresholds(ii)); % this is the hotspot
        OnTarget(1,i,ii) = (sum(VertexSurfaceArea.data(HotSpot&TargetNetwork.data(1:59412,1)==1)) / sum(VertexSurfaceArea.data(HotSpot))); 
        Penalty(1,i,ii) = (sum(VertexSurfaceArea.data(HotSpot&AvoidanceRegion.data(1:59412,1)==1)) / sum(VertexSurfaceArea.data(HotSpot))); 
        
    end
    
end

% clear
% variable
clear MagnE

% save some variables;
save([OutDir '/Optimize/CoilOrientation_OnTarget'],'OnTarget');
save([OutDir '/Optimize/CoilOrientation_Penalty'],'Penalty');

% average accross the e-field thresholds;
AvgOnTarget = mean(OnTarget,3); % on-target;
AvgPenalty = mean(Penalty,3); % penalty;
AvgPenalizedOnTarget = mean(OnTarget,3) - mean(Penalty,3); % on-target - penalty; used for optimization

% read in the template metric file;
G = gifti([OutDir '/Optimize/CoilCenter_OnTarget.shape.gii']);
G.cdata = zeros(size(G.cdata)); % blank slate;

% on-target value;
G.cdata(RefDirsVertices) = AvgOnTarget; % 
save(G,[OutDir '/Optimize/CoilOrientation_OnTarget.shape.gii']); % write out the on-target metric file;
system(['wb_command -metric-dilate ' OutDir '/Optimize/CoilOrientation_OnTarget.shape.gii ' SkinFile ' 2  ' OutDir '/Optimize/CoilOrientation_OnTarget.shape.gii -nearest']);
system(['wb_command -metric-smoothing ' SkinFile ' ' OutDir '/Optimize/CoilOrientation_OnTarget.shape.gii ' num2str(SmoothingFactor) ' ' OutDir '/Optimize/CoilOrientation_OnTarget_s' num2str(SmoothingFactor) '.shape.gii -fix-zeros']);

% penalty value;
G.cdata(RefDirsVertices) = AvgPenalty; % 
save(G,[OutDir '/Optimize/CoilOrientation_Penalty.shape.gii']); % write out the on-target metric file;
system(['wb_command -metric-dilate ' OutDir '/Optimize/CoilOrientation_Penalty.shape.gii ' SkinFile ' 2  ' OutDir '/Optimize/CoilOrientation_Penalty.shape.gii -nearest']);
system(['wb_command -metric-smoothing ' SkinFile ' ' OutDir '/Optimize/CoilOrientation_Penalty.shape.gii ' num2str(SmoothingFactor) ' ' OutDir '/Optimize/CoilOrientation_Penalty_s' num2str(SmoothingFactor) '.shape.gii -fix-zeros']);

% penalized on-target value;
G.cdata(RefDirsVertices) = AvgPenalizedOnTarget; % 
save(G,[OutDir '/Optimize/CoilOrientation_PenalizedOnTarget.shape.gii']); % write out the on-target metric file;
system(['wb_command -metric-dilate ' OutDir '/Optimize/CoilOrientation_PenalizedOnTarget.shape.gii ' SkinFile ' 2  ' OutDir '/Optimize/CoilOrientation_PenalizedOnTarget.shape.gii -nearest']);
system(['wb_command -metric-smoothing ' SkinFile ' ' OutDir '/Optimize/CoilOrientation_PenalizedOnTarget.shape.gii ' num2str(SmoothingFactor) ' ' OutDir '/Optimize/CoilOrientation_PenalizedOnTarget_s' num2str(SmoothingFactor) '.shape.gii -fix-zeros']);

% find the coil orientation that has
% the relative on-target score (on average);
Idx = find(AvgPenalizedOnTarget==max(AvgPenalizedOnTarget));
OrientationVertex = RefDirsVertices(Idx(1)); % this is the final coil placement site;
CoilOrientationCoords = Skin.vertices(OrientationVertex,:);
system(['echo ' num2str(CoilOrientationCoords(1)) ' ' num2str(CoilOrientationCoords(2)) ' ' num2str(CoilOrientationCoords(3)) ' > ' OutDir '/Optimize/CoilOrientationCoordinates.txt']);

% make a foci on the skin surface;
system(['echo Target > ' OutDir '/Optimize/tmp.txt']);
system(['echo 1 0 0 ' num2str(round(CoilOrientationCoords(1))) ' ' num2str(round(CoilOrientationCoords(2))) ' ' num2str(round(CoilOrientationCoords(3))) ' >> ' OutDir '/Optimize/tmp.txt']);
system(['wb_command -foci-create  ' OutDir '/Optimize/CoilOrientation.foci -class CoilOrientation ' OutDir '/Optimize/tmp.txt ' SkinFile]);
system(['rm ' OutDir '/Optimize/tmp.txt']); % remove intermediate file;

% write out the final E-field strength map
MagnE = ft_read_cifti_mod([OutDir '/Optimize/magnE_BestCoilCenter+AllOrientations.dtseries.nii']); MagnE.data = MagnE.data(:,Idx(1));
ft_write_cifti_mod([OutDir '/Optimize/magnE_BestCoilCenter+BestOrientation.dtseries.nii'],MagnE);

% define the path and name for the simulation .log
Log = dir([OutDir '/Optimize/Simulation/*.log']);

% extract the affine transforms already created by SimNIBS;
M = ExtractMatsFromLogFile([Log.folder '/' Log.name]);

% read in the .stl version of the coil used;
Coil = stlread(strrep(CoilModel,'.ccd','.stl'));

% infer the coil name;
CoilName = strsplit(CoilModel,{'/','.'});
CoilName = CoilName{end-1};

% write out the affine transformation matrix;
writematrix(M{Idx(1)},[OutDir '/Optimize/'...
CoilName '_xfm.txt'],'Delimiter','tab');

% convert to .surf.gii
G = gifti; % preallocate;
G.mat = eye(4); % identity matrix
G.vertices = single(Coil.Points); % vertices
G.faces = int32(Coil.ConnectivityList); % edges
save(G,[OutDir '/Optimize/' CoilName '.surf.gii']);

% apply transformation matrix and set structure type;
system(['wb_command -surface-apply-affine ' OutDir '/Optimize/' CoilName '.surf.gii ' OutDir '/Optimize/' CoilName '_xfm.txt ' OutDir '/Optimize/' CoilName '.surf.gii']);
system(['wb_command -set-structure ' OutDir '/Optimize/' CoilName '.surf.gii CEREBELLUM -surface-type RECONSTRUCTION']); % note: calling the TMS coil cerebellum is an ugly workaround, but there is not TMS coil option, of course.

% remove large  intermediate files 
system(['rm -rf ' OutDir '/Optimize/Simulation/']);

end

% TANS sub-functions
function [Output] = SortCircleCoords(Coords,W)
% cjl;

% preallocate;
SortedCoords = zeros(size(Coords));

% preallocate;
Used = [];
Idx = 1;

% sweep the coordinates;
for i = 1:size(Coords,1)
    D = pdist2(Coords(Idx,:),Coords);
    D(Used)=nan;
    Idx = find(D==min(nonzeros(D)));
    Idx = Idx(1); % in case there is more than one;
    SortedCoords(i,:) = Coords(Idx,:);
    Used = [Used Idx];
end

% preallocate;
Output = zeros(size(Coords));

% sweep the coordinates
for i = 0:size(SortedCoords,1)-1
    Idx = circshift(1:size(SortedCoords,1),W+i);
    Output(i+1,:) = mean(SortedCoords(Idx(1:W*2),:));
end

end
function [Output] = ExtractMatsFromLogFile(File)

% setup the import Options
opts = delimitedTextImportOptions("NumVariables", 1);

% specify range and delimiter
opts.DataLines = [1, Inf];
opts.Delimiter = "";

% specify column names and types
opts.VariableTypes = "char";
opts = setvaropts(opts, 1, "WhitespaceRule", "preserve");
opts = setvaropts(opts, 1, "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
Log = readtable(File, opts);

% preallocate
count = 0;

% sweep the lines
% of the log file;
for i = 1:size(Log,1)
    
    tmp = table2cell(Log(i,:));
    tmp = tmp{1};
    
    % if this line represents the first 
    % line of an affine transformation matrix 
    if ~isempty(tmp) && strcmp(tmp(1:2),'[[') 
        
        Idx = i:i+3; % these rows correspond to matrix i
        tmp = table2array(Log(Idx,:));
        tmp = strrep(tmp,'[','');
        tmp = strrep(tmp,']','');
        
        % preallocate;
        m = zeros(4);
        
        % sweep
        % the rows;
        for ii = 1:4
            m(ii,:) = str2num(tmp{ii,1});
        end
        
        % log the matrix
        count = count + 1;
        Output{count} = m;
        
    end
    
end

end




