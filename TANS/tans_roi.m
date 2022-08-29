function tans_roi(TargetNetwork,MidthickSurfs,VertexSurfaceArea,Sulc,SearchSpace,OutDir,Paths)
% cjl; cjl2007@med.cornell.edu;

% Description of input parameters
% "TargetNetwork": Cifti file containing only ttarget network vertices).
% "MidthickSurfs":  Paths to low (32k) dimensional FS_LR midthickness surfaces (Cell array of strings, MidthickSurfs{1} = path to LH, MidthickSurfs{2} = path to RH).
% "VertexSurfaceArea": A CIFTI file containing the vertex surface areas (structure array). 
% "Sulc": CIFTI file containing sulcal depth information (structure array). 
% "SearchSpace": A CIFTI file containing binarized mask used to define the search space. 
% "OutDir": Path to the output folder (string).
% "Paths": Paths to folders that must be added to Matlab search path (cell array of strings). 

% add some directories 
% to the search path;
for i = 1:length(Paths)
addpath(genpath(Paths{i})); % 
end

rng(44); % for reproducibility;

% load midthickness surfaces
LH = gifti(MidthickSurfs{1});
RH = gifti(MidthickSurfs{2});

% extract coordinates for all cortical vertices
BrainStructure = TargetNetwork.brainstructure; % barrow the brain structure index
SurfaceCoordinates = [LH.vertices; RH.vertices]; % combine hemipsheres
SurfaceIndex = BrainStructure > 0 & BrainStructure < 3; % cortex index
BrainStructure(BrainStructure==-1) = []; % remove medial wall vertices
SurfaceIndex = SurfaceIndex(1:size(SurfaceCoordinates,1));
SurfaceCoordinates = SurfaceCoordinates(SurfaceIndex,:);

% find vertices on the medial surface of the brain
D = pdist2(SurfaceCoordinates,SurfaceCoordinates);
D(BrainStructure==1,BrainStructure==1)=nan;
D(BrainStructure==2,BrainStructure==2)=nan;
MedialWallVertices = find(min(D,[],2) < 10); % 10mm seems to work okay; the idea is to this needed to be improved in future

% make the ROI dir.;
mkdir([OutDir '/ROI']);

% write out the search space (useful for double 
% checking that you specified the cortical zones you intended).
ft_write_cifti_mod([OutDir '/ROI/SearchSpace'],SearchSpace);

% this is the full network target (before any editing);
TargetNetwork.data = TargetNetwork.data~=0; % note: should be a binarized map already, but just in case it is not. 

O = TargetNetwork; % preallocate
O.data = zeros(size(TargetNetwork.data)); % blank slate 
O.data(TargetNetwork.data==1) = 1;  % log network 
O.data(59413:end) = 0; % cortex only
ft_write_cifti_mod([OutDir '/ROI/TargetNetwork'],O); % write out the .dtseries.nii;

% generate borders file for target network.
system(['echo Target Network > ' OutDir '/ROI/Labels.txt']);
system(['echo 1 0 0 0 255 >> ' OutDir '/ROI/Labels.txt']);
system(['wb_command -cifti-label-import ' OutDir '/ROI/TargetNetwork.dtseries.nii ' OutDir '/ROI/Labels.txt ' OutDir '/ROI/TargetNetwork.dlabel.nii -discard-others']);
system(['wb_command -cifti-label-to-border ' OutDir '/ROI/TargetNetwork.dlabel.nii -border ' MidthickSurfs{1} ' ' OutDir '/ROI/TargetNetwork.L.border']);
system(['wb_command -cifti-label-to-border ' OutDir '/ROI/TargetNetwork.dlabel.nii -border ' MidthickSurfs{2} ' ' OutDir '/ROI/TargetNetwork.R.border']);
system(['rm ' OutDir '/ROI/Labels.txt ' OutDir '/ROI/TargetNetwork.dlabel.nii']); % remove intermediate files

% this is the full network target (constrained to search space);
TargetNetwork.data(SearchSpace.data==0) = 0; % constrain to search space
O = TargetNetwork; % preallocate
O.data = zeros(size(TargetNetwork.data)); % blank slate 
O.data(TargetNetwork.data==1) = 1;  % log network within the search space 
O.data(59413:end) = 0; % cortex only
ft_write_cifti_mod([OutDir '/ROI/TargetNetwork+SearchSpace'],O); % write out the .dtseries.nii;

% this is the full network target,constrained to search space, after sulcal + medial wall masking;
TargetNetwork.data(Sulc.data < 0) = 0; % remove network vertices in sulcus / fundus;
TargetNetwork.data(MedialWallVertices) = 0; % remove medial surface vertices
O = TargetNetwork; % preallocate
O.data = zeros(size(TargetNetwork.data)); % blank slate 
O.data(TargetNetwork.data==1) = 1;  % log network within the search space 
O.data(59413:end) = 0; % cortex only
ft_write_cifti_mod([OutDir '/ROI/TargetNetwork+SearchSpace+SulcalMask'],O); % write out the .dtseries.nii;

% evaluate the size of each cluster;
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
        ClusterSize(i) =  sum(VertexSurfaceArea.data(find(Clusters.data==i)));
        
    end
    
    % find the largest cluster;
    TargetPatch = uClusters(ClusterSize==max(ClusterSize)); % bigger is better;
    TargetPatch = Clusters.data(1:59412)==TargetPatch;
    
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
    TargetPatch = uClusters(1); % 
    TargetPatch = Clusters.data(1:59412)==TargetPatch;
    
end

% write out the target network patch 
O.data = zeros(size(O.data,1),1); % blank slate
O.data(TargetPatch==1)=1; % binarize;
ft_write_cifti_mod([OutDir '/ROI/TargetNetworkPatch'],O);

end

