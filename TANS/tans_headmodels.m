function tans_headmodels(Subdir,OutDir,Paths)
% cjl; cjl2007@med.cornell.edu;
%
% Inputs
% "Subdir": The path to the subject's folder.
% "OutDir": The path to the output folder.

% define some directories;
addpath(genpath(Paths{1})); % define the path to SimNibs software
addpath(genpath(Paths{2})); % define the path to the folder containing "ft_read_cifti" / "gifti" functions

rng(44); % for reproducibility;

% here are some 
% (kind of) arbitrary
% but hard set values ;
SmoothingFactor = 0.25;
nIterations = 500;

% infer subject name;
str = strsplit(Subdir,'/');
Subject = str{end};
 
% directory;
mkdir(OutDir);

% make a subject folder;
mkdir([OutDir '/HeadModel/']);
cd([OutDir '/HeadModel/']);

% copy T1w image;
system(['cp ' Subdir '/anat/T1w/T1w_acpc_dc_restore.nii.gz '...
OutDir '/HeadModel/T1w_acpc_dc_restore.nii.gz']); 

% check to see if the T2w image exists and proceed appropriately;
if exist([Subdir '/anat/T1w/T2w_acpc_dc_restore.nii.gz'],'file')
    
    % copy the T2w image;
    system(['cp ' Subdir '/anat/T1w/T2w_acpc_dc_restore.nii.gz '...
    OutDir '/HeadModel/T2w_acpc_dc_restore.nii.gz']);
    
    % construct a tetrahedral headmesh using headreco using both the T1w and T2w image;
    system(['headreco all ' Subject ' T1w_acpc_dc_restore.nii.gz T2w_acpc_dc_restore.nii.gz --cat_print --skip_coreg > /dev/null 2>&1']);
    
else
    
    % construct a tetrahedral headmesh using headreco using only the T1w image;
    system(['headreco all ' Subject ' T1w_acpc_dc_restore.nii.gz --cat_print --skip_coreg > /dev/null 2>&1']);
    
end
 
% convert the .stl file to .surf.gii
system(['mris_convert ' OutDir '/HeadModel/m2m_' Subject '/skin.stl '...
OutDir '/HeadModel/m2m_' Subject '/Skin.surf.gii > /dev/null 2>&1']);

% calculate the surface area that each skin mesh vertex is responsible for;
system(['wb_command -surface-vertex-areas ' OutDir '/HeadModel/m2m_' Subject '/Skin.surf.gii '...
OutDir '/HeadModel/m2m_' Subject '/Skin.va.shape.gii']);

% apply some spatial smoothing;
system(['wb_command -surface-smoothing ' OutDir '/HeadModel/m2m_' Subject '/Skin.surf.gii ' num2str(SmoothingFactor) ' ' num2str(nIterations) ' ' OutDir '/HeadModel/m2m_' Subject '/SkinSmoothed.surf.gii']);
system(['wb_command -set-structure ' OutDir '/HeadModel/m2m_' Subject '/SkinSmoothed.surf.gii CORTEX_LEFT']);
system(['wb_command -set-structure ' OutDir '/HeadModel/m2m_' Subject '/Skin.surf.gii CORTEX_LEFT']);

% now, create a distance matrix;

% load the original (non-smoothed) skin surface;
Skin = gifti([OutDir '/HeadModel/m2m_' Subject '/Skin.surf.gii']);

% preallocate the distance matrix;
D = zeros(size(Skin.vertices,1));

% sweep through all the vertices;
for i = 1:size(Skin.vertices,1)
    
    % calculate geodesic distances from vertex i to all other vertices;
    system(['wb_command -surface-geodesic-distance ' OutDir '/HeadModel/m2m_' Subject '/Skin.surf.gii ' num2str(i-1) ' ' OutDir '/HeadModel/tmp_' num2str(i) '.shape.gii']);
    tmp = gifti([OutDir '/HeadModel/tmp_' num2str(i) '.shape.gii']);
    system(['rm ' OutDir '/HeadModel/tmp_' num2str(i) '.shape.gii']);
    D(:,i) = tmp.cdata; % log distances
        
end

% save distance matrix;
save([OutDir '/HeadModel/m2m_' Subject...
'/SkinDistanceMatrix'],'D','-v7.3');

end


