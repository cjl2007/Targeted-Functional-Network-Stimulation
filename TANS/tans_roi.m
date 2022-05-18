function [ROI,LargestCluster] = tans_roi(Subdir,N,SearchSpace,FS,OutDir,Paths)
% cjl; cjl2007@med.cornell.edu;
%
% Inputs
% "Subdir" (Subject Directory): Path to subject folder. 

% define some directories;
addpath(genpath(Paths{1})); % define the path to SimNibs software
addpath(genpath(Paths{2})); % define the path to the folder containing "ft_read_cifti" / "gifti" functions

rng(44); % for reproducibility;

% infer subject name;
str = strsplit(Subdir,'/');
Subject = str{end};

% load FreeSurfer (FS) gyral labels 
% and sulcal depth (Sulc) information;
FS = ft_read_cifti_mod(FS);
Sulc = ft_read_cifti_mod([Subdir '/anat/MNINonLinear/fsaverage_LR32k/' Subject '.sulc.32k_fs_LR.dscalar.nii']);
BrainStructure = FS.brainstructure; % barrow the brain structure index
Sulc.data(BrainStructure==-1)=[]; % remove medial wall vertices
BrainStructure(BrainStructure==-1)=[]; % remove medial wall vertices

% define midthickness surfaces and vertex surface area information;
MidthickSurfs{1} = [Subdir '/anat/T1w/fsaverage_LR32k/' Subject '.L.midthickness.32k_fs_LR.surf.gii'];
MidthickSurfs{2} = [Subdir '/anat/T1w/fsaverage_LR32k/' Subject '.R.midthickness.32k_fs_LR.surf.gii'];
VA = ft_read_cifti_mod([Subdir '/anat/T1w/fsaverage_LR32k/' Subject '.midthickness_va.32k_fs_LR.dscalar.nii']);

% load midthickness surfaces
LH = gifti(MidthickSurfs{1});
RH = gifti(MidthickSurfs{2});

% extract coordinates for all cortical vertices
SurfaceCoordinates = [LH.vertices; RH.vertices]; % combine hemipsheres
surf_indices_incifti = N.brainstructure > 0 & N.brainstructure < 3;
surf_indices_incifti = surf_indices_incifti(1:size(SurfaceCoordinates,1));
SurfaceCoordinates = SurfaceCoordinates(surf_indices_incifti,:);

% find vertices on the medial surface of the brain
D = pdist2(SurfaceCoordinates,SurfaceCoordinates);
D(BrainStructure==1,BrainStructure==1)=nan;
D(BrainStructure==2,BrainStructure==2)=nan;
MedialWallVertices = find(min(D,[],2) < 10); % 10mm seems to work okay;

% make the ROI dir.;
mkdir([OutDir '/ROI']);

% this is the full network 
% target (before any editing);
TargetNetwork = N.data~=0;

O = N; % preallocate
O.data = zeros(size(N.data)); % blank slate 
O.data(TargetNetwork==1) = 1;  % log network 
O.data(59413:end) = 0; % cortex only
ft_write_cifti_mod([OutDir '/ROI/TargetNetwork'],O); % write out the .dtseries.nii;
ROI = O.data(1:59412); % this is the ROI we will pass to subsequent functions;

% write out the first community.
system(['echo Target Network > ' OutDir '/ROI/Labels.txt']);
system(['echo 1 0 0 0 255 >> ' OutDir '/ROI/Labels.txt']);
system(['wb_command -cifti-label-import ' OutDir '/ROI/TargetNetwork.dtseries.nii ' OutDir '/ROI/Labels.txt ' OutDir '/ROI/TargetNetwork.dlabel.nii -discard-others']);
system(['wb_command -cifti-label-to-border ' OutDir '/ROI/TargetNetwork.dlabel.nii -border ' MidthickSurfs{1} ' ' OutDir '/ROI/TargetNetwork.L.border']);
system(['wb_command -cifti-label-to-border ' OutDir '/ROI/TargetNetwork.dlabel.nii -border ' MidthickSurfs{2} ' ' OutDir '/ROI/TargetNetwork.R.border']);
system(['rm ' OutDir '/ROI/Labels.txt ' OutDir '/ROI/TargetNetwork.dlabel.nii']); % remove intermediate files

% this is the full network target (constrained to search space);
TargetNetwork(~ismember(FS.data,SearchSpace)) = 0; % constrain to search space
O = N; % preallocate
O.data = zeros(size(N.data)); % blank slate 
O.data(TargetNetwork==1) = 1;  % log network within the search space 
O.data(59413:end) = 0; % cortex only
ft_write_cifti_mod([OutDir '/ROI/TargetNetwork+SearchSpace'],O); % write out the .dtseries.nii;

% this is the full network target,constrained to search space, after sulcal/medial wall masking;
TargetNetwork(Sulc.data < 0) = 0; % remove network vertices in sulcus / fundus;
TargetNetwork(MedialWallVertices) = 0; % remove medial surface vertices
O = N; % preallocate
O.data = zeros(size(N.data)); % blank slate 
O.data(TargetNetwork==1) = 1;  % log network within the search space 
O.data(59413:end) = 0; % cortex only
ft_write_cifti_mod([OutDir '/ROI/TargetNetwork+SearchSpace+SulcalMask'],O); % write out the .dtseries.nii;

% find clusters;
system(['wb_command -cifti-find-clusters ' OutDir '/ROI/TargetNetwork+SearchSpace+SulcalMask.dtseries.nii 0 0 0 0 COLUMN ' OutDir...
'/ROI/TargetNetwork+SearchSpace+SulcalMask+Clusters.dtseries.nii -left-surface ' MidthickSurfs{1} ' -right-surface ' MidthickSurfs{2}]);

% read in the cluster metric file;
Clusters = ft_read_cifti_mod([OutDir '/ROI/TargetNetwork+SearchSpace+SulcalMask+Clusters.dtseries.nii']); 
uClusters = unique(nonzeros(Clusters.data)); % unique clusters
nClusters = length(uClusters); % count the number of clusters

% if there is more
% than one cluster;
if nClusters > 1
    
    % preallocate;
    ClusterSize = zeros(nClusters,1);
    
    % sweep the clusters;
    for i = 1:nClusters
        
        % calculate the total surface area of this cluster;
        ClusterSize(i) =  sum(VA.data(find(Clusters.data==i)));
        
    end
    
    % find the largest cluster;
    LargestCluster = uClusters(ClusterSize==max(ClusterSize)); % bigger is better;
    LargestCluster = Clusters.data(1:59412)==LargestCluster;
    
    % preallocate 
    % the output;
    O = Clusters;
    O.data = zeros(size(O.data));
    
    % sweep the clusters;
    for i = 1:nClusters
    O.data(Clusters.data==i,1)   = ClusterSize(i); 
    end
    
    % write out the cluster sizes;
    ft_write_cifti_mod([OutDir...
    '/ROI/TargetNetwork+SearchSpace+SulcalMask+Clusters'],O);
    
    
else
    
    % find the largest cluster;
    LargestCluster = uClusters(1); % 
    LargestCluster = Clusters.data(1:59412)==LargestCluster;
    
end

end

