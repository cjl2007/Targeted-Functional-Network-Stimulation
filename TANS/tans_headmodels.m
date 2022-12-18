function tans_headmodels(Subject,T1w,T2w,OutDir,Paths)
% cjl; cjl2007@med.cornell.edu;

% Description of input parameters
% "Subject": The subject ID (string).
% "T1w": Path to T1w-weighted anatomical image (string).
% "T2w": Path to T2w-weighted anatomical image (string).
% "OutDir": Path to the output folder (string).
% "Paths": Paths to folders that must be added to Matlab search path (cell array of strings).  

% Notes: Head model construction benefits from, but does not require, a T2-weighted anatomical image. 
% If a T2-weighted image is not available, T2w can be set to []. If a T2w image is provided, 
% it should already be co-registered to the T1w. 

% Notes: apply a small amount of spatial smoothing; 
SmoothingStrength = 0.50; % Amount of spatial smoothing to apply to skin mesh [0.0 - 1.0].
SmoothingIterations = 10; % The number of iterations when smoothing the skin mesh.

% add some directories 
% to the search path;
for i = 1:length(Paths)
addpath(genpath(Paths{i})); % 
end

% make a subject folder;
mkdir([OutDir '/HeadModel/']);
cd([OutDir '/HeadModel/']);

% if T2w image
% is available;
if ~isempty(T2w)
    
    % copy the T1w & T2w images to Head Model directory;
    system(['cp ' T1w ' ' OutDir '/HeadModel/T1w.nii.gz']); 
    system(['cp ' T2w ' ' OutDir '/HeadModel/T2w.nii.gz']); 

    % construct a tetrahedral headmesh using headreco using both the T1w and T2w images;
    system(['charm ' Subject ' T1w.nii.gz T2w.nii.gz --skipregisterT2 > /dev/null 2>&1']);
    
else
    
    % copy the T1w image to Head Model directory;
    system(['cp ' T1w ' ' OutDir '/HeadModel/T1w.nii.gz']); 

    % construct a tetrahedral headmesh using headreco using only the T1w image;
    system(['charm ' Subject ' T1w.nii.gz --skipregisterT2 > /dev/null 2>&1']);
    
end

% extract the skin tissue and "background" compartment; 
system(['mri_binarize --i ' OutDir '/HeadModel/m2m_' Subject '/final_tissues.nii.gz --match 5 --o ' OutDir '/HeadModel/m2m_' Subject '/skin.nii.gz > /dev/null 2>&1']);
system(['mri_binarize --i ' OutDir '/HeadModel/m2m_' Subject '/segmentation/labeling.nii.gz --match 517 0 --o ' OutDir '/HeadModel/m2m_' Subject '/background.nii.gz --inv > /dev/null 2>&1 ']);
system(['fslmaths ' OutDir '/HeadModel/m2m_' Subject '/skin.nii.gz -add ' OutDir '/HeadModel/m2m_' Subject '/background.nii.gz -bin ' OutDir '/HeadModel/m2m_' Subject '/skin.nii.gz > /dev/null 2>&1']);

% remove "edge" effects 
% (otherwise we can end up with holes in our skin mesh)
nii = niftiread([OutDir '/HeadModel/m2m_' Subject '/skin.nii.gz']);
dims = size(nii); % nii dimensions
nii([1 dims(1)],:,:) = 0; 
nii(:,[1 dims(2)],:) = 0; 
nii(:,:,[1 dims(3)]) = 0;
nii_info = niftiinfo([OutDir '/HeadModel/m2m_' Subject '/skin.nii.gz']);
niftiwrite(nii,[OutDir '/HeadModel/m2m_' Subject '/skin_noedges'],nii_info);
system(['gzip ' OutDir '/HeadModel/m2m_' Subject '/skin_noedges.nii -f']);

% create skin surface mesh;
system(['mri_tessellate -n ' OutDir '/HeadModel/m2m_' Subject '/skin_noedges.nii.gz 1 ' OutDir '/HeadModel/m2m_' Subject '/skin.orig > /dev/null 2>&1']);
system(['mris_convert ' OutDir '/HeadModel/m2m_' Subject '/skin.orig ' OutDir '/HeadModel/m2m_' Subject '/Skin.surf.gii > /dev/null 2>&1']);
system(['wb_command -surface-smoothing ' OutDir '/HeadModel/m2m_' Subject '/Skin.surf.gii ' num2str(SmoothingStrength) ' ' num2str(SmoothingIterations) ' ' OutDir '/HeadModel/m2m_' Subject '/Skin.surf.gii > /dev/null 2>&1']);
system(['wb_command -set-structure ' OutDir '/HeadModel/m2m_' Subject '/Skin.surf.gii CORTEX_LEFT -surface-type RECONSTRUCTION > /dev/null 2>&1']);

% now, prepare some files that can
% be imported into BrainSight later on if desired,

% create a .stl file as well
system(['mris_convert ' OutDir '/HeadModel/m2m_' Subject '/Skin.surf.gii '...
OutDir '/HeadModel/m2m_' Subject '/Skin.stl > /dev/null 2>&1']);

% load the head model for this subject;
M = mesh_load_gmsh4([OutDir '/HeadModel/m2m_' Subject '/' Subject '.msh']);
S = mesh_extract_regions(M,'region_idx',[2 1002]); % 2 == GM.

% convert to .surf.gii
G = gifti; % preallocate;
G.mat = eye(4); % identity matrix
G.vertices = single(S.nodes); % vertices
G.faces = int32(S.triangles); % edges
save(G,[OutDir '/HeadModel/m2m_' Subject '/GrayMatter.surf.gii']);

% apply spatial smoothing and write out as a .stl file;
system(['wb_command -surface-smoothing ' OutDir '/HeadModel/m2m_' Subject '/GrayMatter.surf.gii ' num2str(SmoothingStrength) ' ' num2str(SmoothingIterations) ' ' OutDir '/HeadModel/m2m_' Subject '/GrayMatter.surf.gii > /dev/null 2>&1']);
system(['mris_convert ' OutDir '/HeadModel/m2m_' Subject '/GrayMatter.surf.gii ' OutDir '/HeadModel/m2m_' Subject '/GrayMatter.stl']);

end


