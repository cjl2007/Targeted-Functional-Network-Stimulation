function [SubSampledSearchGrid] = tans_searchgrid(Subdir,ROI,HeadModelDir,OutDir,GridSpacing,CoiltoCortexDistance,Paths)
% cjl; cjl2007@med.cornell.edu;
%
% Inputs
% "Subdir" (Subject Directory): Path to subject folder. 
% "HeadModelDir" (Head Model Directory): Must match
% "OutDir": Output directory containing initial simulations. 
% "GridSpacing": average distance between points in search grid (in mm).
% Other notes: the script assumes HCP-like directory structures. 

% define some directories;
addpath(genpath(Paths{1})); % define the path to SimNibs software
addpath(genpath(Paths{2})); % define the path to the folder containing "ft_read_cifti" / "gifti" functions

rng(44); % for reproducibility;

% infer subject name;
str = strsplit(Subdir,'/');
Subject = str{end};

% make search grid directory;
mkdir([OutDir '/SearchGrid/']);

% load pial surfaces; 
LH = gifti([Subdir '/anat/T1w/fsaverage_LR32k/' Subject '.L.pial.32k_fs_LR.surf.gii']);
RH = gifti([Subdir '/anat/T1w/fsaverage_LR32k/' Subject '.R.pial.32k_fs_LR.surf.gii']);

% extract coordinates for all cortical vertices 
FS = ft_read_cifti_mod([Subdir '/anat/MNINonLinear/fsaverage_LR32k/' Subject '.aparc.32k_fs_LR.dlabel.nii']); % in principle this could be any cifti file from this subject;
SurfaceCoords = [LH.vertices; RH.vertices]; % combine hemipsheres 
NotMedialWall = FS.brainstructure > 0 & FS.brainstructure < 3;
NotMedialWall = NotMedialWall(1:size(SurfaceCoords,1));
SurfaceCoords = SurfaceCoords(NotMedialWall,:);

% create a search grid above the ROI centroid;

% read in the skin mesh file;
AllSkin = gifti([HeadModelDir '/m2m_'...
Subject '/Skin.surf.gii']);

% identify the ROI centroid, calculate distance to scalp, write out metric file;
system(['wb_command -volume-to-surface-mapping ' HeadModelDir '/m2m_' Subject '/skin.nii.gz ' HeadModelDir '/m2m_' Subject '/Skin.surf.gii ' OutDir '/SearchGrid/DistanceToROI.shape.gii -trilinear']);
G = gifti([OutDir '/SearchGrid/DistanceToROI.shape.gii']); G.cdata = zeros(size(G.cdata)); % blank slate;
G.cdata = pdist2(AllSkin.vertices,mean(SurfaceCoords(ROI==1,:))); % log the distances
save(G,[OutDir '/SearchGrid/DistanceToROI.shape.gii']); % write out the metric file

% generate foci representing scalp 
% vertex directly above target (for visualization purposes);
system(['echo DirectlyAboveTargetVertex > ' OutDir '/SearchGrid/tmp.txt']);
DirectlyAboveTargetCoords = AllSkin.vertices(G.cdata==min(G.cdata),:); % 3-D coordinates
system(['echo 1 0 0 ' num2str(DirectlyAboveTargetCoords(1)) ' ' num2str(DirectlyAboveTargetCoords(2)) ' ' num2str(DirectlyAboveTargetCoords(3)) ' >> ' OutDir '/SearchGrid/tmp.txt']);
system(['wb_command -foci-create  ' OutDir '/SearchGrid/DirectlyAboveTarget.foci -class CoilCenter ' OutDir '/SearchGrid/tmp.txt ' HeadModelDir '/m2m_' Subject '/Skin.surf.gii']);
system(['rm ' OutDir '/SearchGrid/tmp.txt']); % remove intermediate file;

% create a metric file showing the full search grid;
G.cdata(G.cdata>=CoiltoCortexDistance) = 0; G.cdata(G.cdata~=0) = 1;
save(G,[OutDir '/SearchGrid/FullSearchGrid.shape.gii']); % 

% full search grid 
% coordinates & vertices;
SearchGridVertices = find(G.cdata~=0);
FullSearchGrid = AllSkin.vertices(SearchGridVertices,:); % 

% load distance matrix (vertex to vertex distance in geodesic space);
D = smartload([HeadModelDir '/m2m_' Subject '/SkinDistanceMatrix.mat']);

% preallocate neighbors variable;
V = nan(size(SearchGridVertices,1),10);

% sweep through the full search grid
for i = 1:size(SearchGridVertices,1)
    tmp = find(D(SearchGridVertices(i),:) <= GridSpacing);
    V(i,1:length(tmp)) = tmp;
end

% preallocate;
SubSample = [];
Neighbors = [];

% sweep all the
% cortical vertices
for i = 1:size(V,1)
    if ~ismember(SearchGridVertices(i),Neighbors)
        SubSample = [SubSample SearchGridVertices(i)];
        Neighbors = [Neighbors V(i,2:end)];
    end
end

% this is a sensible subsampling of the
% full search grid that is suggested for following analyses
SubSampledSearchGrid = AllSkin.vertices(SubSample,:); 

% create a metric file showing the subsampled search grid
G = gifti([OutDir '/SearchGrid/FullSearchGrid.shape.gii']); 
G.cdata = zeros(size(G.cdata)); % blank slate;
G.cdata(SubSample) = 1:length(SubSample); % log serial numbers
save(G,[OutDir '/SearchGrid/SubSampledSearchGrid.shape.gii']); % write out the metric file

% save some variables;
save([OutDir '/SearchGrid/SubSampledSearchGrid'],'SubSampledSearchGrid');
save([OutDir '/SearchGrid/FullSearchGrid'],'FullSearchGrid');

end