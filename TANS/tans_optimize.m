function tans_optimize(Subject,TargetNetwork,AvoidanceRegion,PercentileThresholds,SkinVol,SkinSurf,SkinDistanceMatrix,VertexSurfaceArea,MidthickSurfs,WhiteSurfs,PialSurfs,MedialWallMask,HeadMesh,PositionUncertainty,CoilModel,AbsoluteThreshold,DiDt,OutDir,Paths)
% cjl; cjl2007@med.cornell.edu;
%
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
% PositionUncertainty: Coil center positioning uncertainty (in mm, numeric). Optimal coil placement is defined as the position that maximizes the on target value on average within this distance. 
% CoilModel: Name of the coil model (string). Must be available in SimNIBS library ( /SimNIBS-3.2/simnibs/ccd-files/). 
% AbsoluteThreshold: This value (in V/m units) represents an assumed neural activation threshold (numeric). Default is 100 V/m. 
% DiDt: A range of stimulation intensities (the speed of variation of the current throughout the stimulating coil, in units of A/us, numeric). Default is DiDt = linspace(1,155,20) * 1e6.
% OutDir: Path to the output folder (string).
% Paths: Paths to folders that must be added to Matlab search path (cell array of strings). 

% add some directories 
% to the search path;
for i = 1:length(Paths)
addpath(genpath(Paths{i})); % 
end

rng(44); % for reproducibility;
warning ('off','all'); % turn off annoying warnings; users could comment this line out if they want. 
SmoothingFactor = 0.85; % 0.85 works well, might need to increase if the search grid is very sparse.

% evaluate E-fields associated with all the coil placements;

% read in some cifti files we will need;
VertexSurfaceArea = ft_read_cifti_mod(VertexSurfaceArea); % vertex surface areas;
SkinSurf2 = gifti(SkinSurf); % note: SkinSurf == path to file and SkinSurf2 = loaded version of file

% load the previously generated metric file;
G_OnTarget = gifti([OutDir '/SearchGrid/SubSampledSearchGrid.shape.gii']); 
SearchGridVertices = find(G_OnTarget.cdata~=0); % 

% count the number of files (note: some attempts will have failed,
% for whatever reason, so we cant assume that file X is simulation X;
files = dir([OutDir '/SearchGrid/normE_*.dtseries.nii']);

% read the first file & count how many orientations were attempted
NormE = ft_read_cifti_mod([OutDir '/SearchGrid/' files(1).name]);
nCols = size(NormE.data,2);

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
    NormE = ft_read_cifti_mod([OutDir '/SearchGrid/' files(i).name]);
    system(['rm ' OutDir '/SearchGrid/' files(i).name]); % remove intermediate files;

    % make sure we have the correct point in the space grid;
    tmp = strsplit(files(i).name,{'_','.'}); % note: sometimes a given point in the search grid will fail;
    idx = str2double(tmp{2});
    
    % sweep the coil orientations;
    for ii = 1:size(NormE.data,2)
        
        % sweep all of the thresholds;
        for iii = 1:length(PercentileThresholds)
            
            HotSpot = NormE.data(1:59412,ii) > prctile(NormE.data(1:59412,ii),PercentileThresholds(iii)); % this is the hotspot
            OnTarget(idx,ii,iii) = (sum(VertexSurfaceArea.data(HotSpot&TargetNetwork.data(1:59412,1)==1)) / sum(VertexSurfaceArea.data(HotSpot))); 
            Penalty(idx,ii,iii) = (sum(VertexSurfaceArea.data(HotSpot&AvoidanceRegion.data(1:59412,1)==1)) / sum(VertexSurfaceArea.data(HotSpot)));
            
        end
        
    end

    % clear 
    % variable
    clear NormE
    
end

% save some variables;
mkdir([OutDir '/Optimize']); % make the dir.;
save([OutDir '/Optimize/CoilCenter_OnTarget'],'OnTarget');
save([OutDir '/Optimize/CoilCenter_Penalty'],'Penalty');

% average accross the e-field thresholds;
AvgOnTarget = mean(OnTarget,3); % on-target;
AvgPenalty = mean(Penalty,3); % penalty;
AvgPenalizedOnTarget = mean(OnTarget,3) - mean(Penalty,3); % on-target - penalty; relative on-target value used for optimization

% create and read in the template metric file;
system(['wb_command -volume-to-surface-mapping ' SkinVol ' ' SkinSurf ' ' OutDir '/Optimize/CoilCenter_OnTarget.shape.gii -trilinear']); % 
G_OnTarget = gifti([OutDir '/Optimize/CoilCenter_OnTarget.shape.gii']);
G_OnTarget.cdata = zeros(size(G_OnTarget.cdata)); % blank slate;

% write out the on-target metric file;
G_OnTarget.cdata(SearchGridVertices) = max(AvgOnTarget,[],2); % average across orientations, for now.
save(G_OnTarget,[OutDir '/Optimize/CoilCenter_OnTarget.shape.gii']); % write out the on-target metric file;
system(['wb_command -metric-smoothing ' SkinSurf ' ' OutDir '/Optimize/CoilCenter_OnTarget.shape.gii ' num2str(SmoothingFactor) ' ' OutDir '/Optimize/CoilCenter_OnTarget_s' num2str(SmoothingFactor) '.shape.gii -fix-zeros']);
G_OnTarget = gifti([OutDir '/Optimize/CoilCenter_OnTarget_s' num2str(SmoothingFactor) '.shape.gii']); % read in the smoothed file;

% write out the penalty metric file;
G_Penalty = G_OnTarget; % preallocate
G_Penalty.cdata = zeros(size(G_OnTarget.cdata)); % blank slate;
G_Penalty.cdata(SearchGridVertices) = max(AvgPenalty,[],2); % average across orientations, for now.
save(G_Penalty,[OutDir '/Optimize/CoilCenter_Penalty.shape.gii']); % write out the on-target metric file;
system(['wb_command -metric-smoothing ' SkinSurf ' ' OutDir '/Optimize/CoilCenter_Penalty.shape.gii ' num2str(SmoothingFactor) ' ' OutDir '/Optimize/CoilCenter_Penalty_s' num2str(SmoothingFactor) '.shape.gii -fix-zeros']);

% write out the penalized on-target metric file;
G_PenalizedOnTarget = G_OnTarget; % preallocate
G_PenalizedOnTarget.cdata = zeros(size(G_OnTarget.cdata)); % blank slate;
G_PenalizedOnTarget.cdata(SearchGridVertices) = max(AvgPenalizedOnTarget,[],2); % average across orientations, for now.
save(G_PenalizedOnTarget,[OutDir '/Optimize/CoilCenter_PenalizedOnTarget.shape.gii']); % write out the relative on-target metric file;
system(['wb_command -metric-smoothing ' SkinSurf ' ' OutDir '/Optimize/CoilCenter_PenalizedOnTarget.shape.gii ' num2str(SmoothingFactor) ' ' OutDir '/Optimize/CoilCenter_PenalizedOnTarget_s' num2str(SmoothingFactor) '.shape.gii -fix-zeros']);
G_PenalizedOnTarget = gifti([OutDir '/Optimize/CoilCenter_PenalizedOnTarget_s' num2str(SmoothingFactor) '.shape.gii']); % read in the smoothed file;

% read in the skin distance matrix;
DistanceMatrix = smartload(SkinDistanceMatrix);

% preallocate the overall quality distance weighted score;
DistanceWeightedScore = zeros(length(SearchGridVertices),1);

% sweep the search grid vertices;
for i = 1:length(SearchGridVertices)
    DistanceWeightedScore(i) = mean(G_PenalizedOnTarget.cdata(DistanceMatrix(SearchGridVertices(i),:) <= PositionUncertainty)); % average value score of all vertices within the specified distance of vertex "i" 
end

% save the the overall quality distance weighted score;
save([OutDir '/Optimize/CoilCenter_DistanceWeightedScore_' num2str(PositionUncertainty) 'mm'],...
'DistanceWeightedScore');

% find the coil placement that has the best quality score (on avg.)
% across vertices within the specified distance (5 mm by default);
Idx = find(DistanceWeightedScore==max(DistanceWeightedScore )); 
CoilCenterVertex = SearchGridVertices(Idx(1)); % this is the final coil placement site;
CoilCenterCoords = SkinSurf2.vertices(CoilCenterVertex,:);
system(['echo ' num2str(CoilCenterCoords(1)) ' ' num2str(CoilCenterCoords(2)) ' ' num2str(CoilCenterCoords(3)) ' > ' OutDir '/Optimize/CoilCenterCoordinates.txt']); % write coordinates out to .txt file;

% make a foci on the skin surface;
system(['echo Target > ' OutDir '/Optimize/tmp.txt']);
system(['echo 1 0 0 ' num2str(CoilCenterCoords(1)) ' ' num2str(CoilCenterCoords(2)) ' ' num2str(CoilCenterCoords(3)) ' >> ' OutDir '/Optimize/tmp.txt']);
system(['wb_command -foci-create  ' OutDir '/Optimize/CoilCenter.foci -class CoilCenter ' OutDir '/Optimize/tmp.txt ' SkinSurf]);
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
D = pdist2(SkinSurf2.vertices,CoilCenterCoords);
A = find(D < 19); % inner diameter;
B = find(D < 20); % outer diameter;
% note: after some trial and error, we found 19/20 works well. This could
% be improved in the future. 

% sort coordinates in a radial fashion (this time with no subsampling);
[RefDirs] = SortCircleCoords(SkinSurf2.vertices(B(~ismember(B,A)),:)); 

% preallocate;
RefDirsVertices = ...
zeros(size(RefDirs,1),1);

% sweep a range of
% coil orientations;
for i = 1:length(RefDirs)
    
    % specify & save coil position;
    s.poslist{1}.pos(i).centre = CoilCenterCoords;
    s.poslist{1}.pos(i).pos_ydir = [RefDirs(i,1), RefDirs(i,2), RefDirs(i,3)];
    
    % while we are at it, lets log 
    % the corresponding skin vertex numbers;
    Tmp = pdist2(SkinSurf2.vertices,RefDirs(i,:));
    Idx = find(Tmp==min(Tmp));
    RefDirsVertices(i) = Idx(1);
    
end
    
% write to volume;
s.map_to_vol = true;
s.fields = 'e'; % normE only;

% remove the dir. if it exists;
if exist([OutDir '/Optimize/Simulation/'],'dir')
system(['rm -rf ' OutDir '/Optimize/Simulation/']);
end

%run simulation;
run_simnibs(s);

% map concatenated volume to the 32k surface;
system(['fslmerge -t ' OutDir '/Optimize/Simulation/subject_volumes/normE.nii.gz ' OutDir '/Optimize/Simulation/subject_volumes/' Subject '*_normE.nii.gz']);
system(['rm ' OutDir '/Optimize/Simulation/subject_volumes/' Subject '*_normE.nii.gz']); % remove intermediate files;
system(['wb_command -volume-to-surface-mapping ' OutDir '/Optimize/Simulation/subject_volumes/normE.nii.gz  ' MidthickSurfs{1} ' ' OutDir '/Optimize/Simulation/subject_volumes/normE.L.32k_fs_LR.shape.gii -ribbon-constrained ' WhiteSurfs{1} ' ' PialSurfs{1} ' -interpolate ENCLOSING_VOXEL']);
system(['wb_command -volume-to-surface-mapping ' OutDir '/Optimize/Simulation/subject_volumes/normE.nii.gz  ' MidthickSurfs{2} ' ' OutDir '/Optimize/Simulation/subject_volumes/normE.R.32k_fs_LR.shape.gii -ribbon-constrained ' WhiteSurfs{2} ' ' PialSurfs{2} ' -interpolate ENCLOSING_VOXEL']);
system(['wb_command -metric-mask ' OutDir '/Optimize/Simulation/subject_volumes/normE.L.32k_fs_LR.shape.gii ' MedialWallMask{1} ' ' OutDir '/Optimize/Simulation/subject_volumes/normE.L.32k_fs_LR.shape.gii']);
system(['wb_command -metric-mask ' OutDir '/Optimize/Simulation/subject_volumes/normE.R.32k_fs_LR.shape.gii ' MedialWallMask{2} ' ' OutDir '/Optimize/Simulation/subject_volumes/normE.R.32k_fs_LR.shape.gii']);
system(['wb_command -cifti-create-dense-timeseries ' OutDir '/Optimize/Simulation/subject_volumes/normE.dtseries.nii -left-metric ' OutDir '/Optimize/Simulation/subject_volumes/normE.L.32k_fs_LR.shape.gii -roi-left ' MedialWallMask{1} ' -right-metric ' OutDir '/Optimize/Simulation/subject_volumes/normE.R.32k_fs_LR.shape.gii -roi-right ' MedialWallMask{2}]);
system(['rm ' OutDir '/Optimize/Simulation/subject_volumes/*shape*']);

% rename the cifti file and remove some intermediate files;
system(['mv ' OutDir '/Optimize/Simulation/subject_volumes/normE.dtseries.nii '...
OutDir '/Optimize/normE_BestCoilCenter+AllOrientations.dtseries.nii']);
system(['rm -rf ' OutDir '/Optimize/Simulation/']);

% read the cifti file & count how many orientations were attempted
NormE = ft_read_cifti_mod([OutDir '/Optimize/normE_BestCoilCenter+AllOrientations.dtseries.nii']);
nCols = size(NormE.data,2);

% preallocate;
OnTarget = zeros(1,nCols,length(PercentileThresholds)); % "OnTaget" variable (% of E-field hotspot that contains target network vertices);
Penalty = zeros(1,nCols,length(PercentileThresholds)); % "Penalty" variable (% of E-field hotspot that contains avoidance region / network vertices);

% sweep the coil orientations;
for i = 1:size(NormE.data,2)
    
    % sweep all of the thresholds;
    for ii = 1:length(PercentileThresholds)

        HotSpot = NormE.data(1:59412,i) > prctile(NormE.data(1:59412,i),PercentileThresholds(ii)); % this is the hotspot
        OnTarget(1,i,ii) = (sum(VertexSurfaceArea.data(HotSpot&TargetNetwork.data(1:59412,1)==1)) / sum(VertexSurfaceArea.data(HotSpot))); 
        Penalty(1,i,ii) = (sum(VertexSurfaceArea.data(HotSpot&AvoidanceRegion.data(1:59412,1)==1)) / sum(VertexSurfaceArea.data(HotSpot))); 
        
    end
    
end

% clear
% variable
clear NormE

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
system(['wb_command -metric-dilate ' OutDir '/Optimize/CoilOrientation_OnTarget.shape.gii ' SkinSurf ' 2  ' OutDir '/Optimize/CoilOrientation_OnTarget.shape.gii -nearest']);
system(['wb_command -metric-smoothing ' SkinSurf ' ' OutDir '/Optimize/CoilOrientation_OnTarget.shape.gii ' num2str(SmoothingFactor) ' ' OutDir '/Optimize/CoilOrientation_OnTarget_s' num2str(SmoothingFactor) '.shape.gii -fix-zeros']);

% penalty value;
G.cdata(RefDirsVertices) = AvgPenalty; % 
save(G,[OutDir '/Optimize/CoilOrientation_Penalty.shape.gii']); % write out the on-target metric file;
system(['wb_command -metric-dilate ' OutDir '/Optimize/CoilOrientation_Penalty.shape.gii ' SkinSurf ' 2  ' OutDir '/Optimize/CoilOrientation_Penalty.shape.gii -nearest']);
system(['wb_command -metric-smoothing ' SkinSurf ' ' OutDir '/Optimize/CoilOrientation_Penalty.shape.gii ' num2str(SmoothingFactor) ' ' OutDir '/Optimize/CoilOrientation_Penalty_s' num2str(SmoothingFactor) '.shape.gii -fix-zeros']);

% penalized on-target value;
G.cdata(RefDirsVertices) = AvgPenalizedOnTarget; % 
save(G,[OutDir '/Optimize/CoilOrientation_PenalizedOnTarget.shape.gii']); % write out the on-target metric file;
system(['wb_command -metric-dilate ' OutDir '/Optimize/CoilOrientation_PenalizedOnTarget.shape.gii ' SkinSurf ' 2  ' OutDir '/Optimize/CoilOrientation_PenalizedOnTarget.shape.gii -nearest']);
system(['wb_command -metric-smoothing ' SkinSurf ' ' OutDir '/Optimize/CoilOrientation_PenalizedOnTarget.shape.gii ' num2str(SmoothingFactor) ' ' OutDir '/Optimize/CoilOrientation_PenalizedOnTarget_s' num2str(SmoothingFactor) '.shape.gii -fix-zeros']);

% find the coil orientation that has
% the relative on-target score (on average);
Idx = find(AvgPenalizedOnTarget==max(AvgPenalizedOnTarget));
OrientationVertex = RefDirsVertices(Idx(1)); % this is the final coil placement site;
CoilOrientationCoords = SkinSurf2.vertices(OrientationVertex,:);
system(['echo ' num2str(CoilOrientationCoords(1)) ' ' num2str(CoilOrientationCoords(2)) ' ' num2str(CoilOrientationCoords(3)) ' > ' OutDir '/Optimize/CoilOrientationCoordinates.txt']);

% make a foci on the skin surface;
system(['echo Target > ' OutDir '/Optimize/tmp.txt']);
system(['echo 1 0 0 ' num2str(CoilOrientationCoords(1)) ' ' num2str(CoilOrientationCoords(2)) ' ' num2str(CoilOrientationCoords(3)) ' >> ' OutDir '/Optimize/tmp.txt']);
system(['wb_command -foci-create  ' OutDir '/Optimize/CoilOrientation.foci -class CoilOrientation ' OutDir '/Optimize/tmp.txt ' SkinSurf]);
system(['rm ' OutDir '/Optimize/tmp.txt']); % remove intermediate file;

% write out the final E-field strength map
NormE = ft_read_cifti_mod([OutDir '/Optimize/normE_BestCoilCenter+AllOrientations.dtseries.nii']); NormE.data = NormE.data(:,Idx(1));
ft_write_cifti_mod([OutDir '/Optimize/normE_BestCoilCenter+BestOrientation.dtseries.nii'],NormE);

%% now, lets estimate the optimal dose (assuming a specific threshold for neural activation).  

% if stimulation intensities 
% to test are not specified  
if isempty(DiDt)
DiDt = linspace(1,155,20) * 1e6; % A/us
end 

% sweep intensity levels
for i = 1:length(DiDt)
    NormE.data(:,i) = (NormE.data(:,1) / (1*1e6)) * DiDt(i); 
end

% write out the cifti file; each column == a different stimulation intensity level
ft_write_cifti_mod([OutDir '/Optimize/normE_BestCoilCenter+BestOrientation_AllStimulationItensities.dtseries.nii'],NormE);
system(['echo ' num2str(DiDt) ' > ' OutDir '/Optimize/DiDt.txt']); % write out the intensities used;

% if no absolute 
% threshold is specified
if isempty(AbsoluteThreshold)
AbsoluteThreshold = 100; % V/m
end

% preallocate;
OnTarget = zeros(size(NormE.data,2),1); % "OnTaget" variable (% of E-field hotspot that contains target network vertices);
Penalty = zeros(size(NormE.data,2),1); % "Penalty" variable (% of E-field hotspot that contains avoidance region / network vertices);

% sweep the range of 
% stimulation intensities
for t = 1:size(NormE.data,2)
    
    HotSpot = NormE.data(1:59412,t) >= AbsoluteThreshold; % this is the hotspot
    
    % only consider if hotspot
    % is at least 1000 mm2
    if sum(VertexSurfaceArea.data(HotSpot)) > 10^3
        
        % calculate proportion of the suprathreshold E-field that is on-target
        OnTarget(t,1) = sum(VertexSurfaceArea.data(HotSpot&TargetNetwork.data(1:59412,1)==1)) / sum(VertexSurfaceArea.data(HotSpot)); 
        Penalty(t,1) = sum(VertexSurfaceArea.data(HotSpot&AvoidanceRegion.data(1:59412,1)==1)) / sum(VertexSurfaceArea.data(HotSpot)); 
 
    end

end

H = figure; % prellocate parent figure
set(H,'position',[1 1 350 225]); hold;

% plot the results;
plot(OnTarget * 100,'Color','k','LineWidth',1); hold;

% make it "pretty";
xlim([0 length(DiDt)]);
ylim([0 100]);
yticks(0:20:100);
xlabel('dI/dt (A/\mus)');
ylabel('% On-Target')
xticks(1:1:length(DiDt));
xticklabels(round(DiDt/1e6));
xtickangle(90);
set(gca,'FontName','Arial','FontSize',12,'TickLength',[0 0]);
saveas(gcf,[OutDir '/Optimize/PercentageOnTarget_AbsoluteThreshold_' num2str(AbsoluteThreshold) '_Vm.pdf']); 
close all;

H = figure; % prellocate parent figure
set(H,'position',[1 1 350 225]); hold;

% plot the results;
plot(Penalty * 100,'Color','k','LineWidth',1); hold;

% make it "pretty";
xlim([0 length(DiDt)]);
ylim([0 100]);
yticks(0:20:100);
xlabel('dI/dt (A/\mus)');
ylabel('Penalty')
xticks(1:1:length(DiDt));
xticklabels(round(DiDt/1e6));
xtickangle(90);
set(gca,'FontName','Arial','FontSize',12,'TickLength',[0 0]);
saveas(gcf,[OutDir '/Optimize/PercentageOffTarget_AbsoluteThreshold_' num2str(AbsoluteThreshold) '_Vm.pdf']); 
close all;

H = figure; % prellocate parent figure
set(H,'position',[1 1 350 225]); hold;

% plot the results;
plot((OnTarget - Penalty) * 100,'Color','k','LineWidth',1); hold;

% make it "pretty";
xlim([0 length(DiDt)]);
ylim([0 100]);
yticks(0:20:100);
xlabel('dI/dt (A/\mus)');
ylabel('% On-Target - Penalty')
xticks(1:1:length(DiDt));
xticklabels(round(DiDt/1e6));
xtickangle(90);
set(gca,'FontName','Arial','FontSize',12,'TickLength',[0 0]);
saveas(gcf,[OutDir '/Optimize/PercentageOnTarget+Penalty_AbsoluteThreshold_' num2str(AbsoluteThreshold) '_Vm.pdf']); 
close all;

end

function [Output] = SortCircleCoords(Coords)
% cjl;

% preallocate;
Output = zeros(size(Coords));

% preallocate;
Used = [];
Idx = 1;

% sweep the coordinates;
for i = 1:size(Coords,1)
    D = pdist2(Coords(Idx,:),Coords);
    D(Used)=nan;
    Idx = find(D==min(nonzeros(D)));
    Idx = Idx(1); % in case there is more than one;
    Output(i,:) = Coords(Idx,:);
    Used = [Used Idx];
end

end

%
function out = smartload(matfile)

out = load(matfile);
names = fieldnames(out);
out = eval(['out.' names{1}]);

end