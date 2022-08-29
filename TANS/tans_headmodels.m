function tans_headmodels(Subject,T1w,T2w,SmoothingStrength,SmoothingIterations,OutDir,Paths)
% cjl; cjl2007@med.cornell.edu;

% Description of input parameters
% "Subject": The subject ID (string).
% "T1w": Path to T1w-weighted anatomical image (string).
% "T2w": Path to T2w-weighted anatomical image (string). 
% "SmoothingStrength": Amount of spatial smoothing to apply to skin mesh [0.0 - 1.0].
% "SmoothingIterations": The number of iterations when smoothing the skin mesh.
% "OutDir": Path to the output folder (string).
% "Paths": Paths to folders that must be added to Matlab search path (cell array of strings).  

% Notes: Headreco benefits from but does not require a T2-weighted anatomical image. 
% If a T2-weighted image is not available, T2w can be set to []. If a T2w image is provided, 
% it should already be co-registered to the T1w. 

% add some directories 
% to the search path;
for i = 1:length(Paths)
addpath(genpath(Paths{i})); % 
end

rng(44); % for reproducibility;

% if no smoothing strength variable 
% is specified, use the default value of 0.25
if isempty(SmoothingStrength)
SmoothingStrength = 0.25; % 
end

% if no smoothing iterations variable 
% is specified, use the default value of 500
if isempty(SmoothingIterations)
SmoothingIterations = 500; % 
end

% Notes: these "default values" for smoothing strength and smoothing
% iterations work well in my data, but not in yours. Ultimately, the
% smoothed surface is usually only used for visualization purposes, so no
% need to worry too much about these settings. 
 
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
    system(['headreco all ' Subject ' T1w.nii.gz T2w.nii.gz --cat_print --skip_coreg > /dev/null 2>&1']);
    
else
    
    % copy the T1w image to Head Model directory;
    system(['cp ' T1w ' ' OutDir '/HeadModel/T1w.nii.gz']); 

    % construct a tetrahedral headmesh using headreco using only the T1w image;
    system(['headreco all ' Subject ' T1w.nii.gz --cat_print --skip_coreg > /dev/null 2>&1']);
    
end
 
% convert the .stl file to .surf.gii
system(['mris_convert ' OutDir '/HeadModel/m2m_' Subject '/skin.stl '...
OutDir '/HeadModel/m2m_' Subject '/Skin.surf.gii > /dev/null 2>&1']);

% calculate the surface area that each vertex is responsible for;
system(['wb_command -surface-vertex-areas ' OutDir '/HeadModel/m2m_' Subject '/Skin.surf.gii '...
OutDir '/HeadModel/m2m_' Subject '/Skin.va.shape.gii']);

% apply some spatial smoothing to skin mesh;
system(['wb_command -surface-smoothing ' OutDir '/HeadModel/m2m_' Subject '/Skin.surf.gii ' num2str(SmoothingStrength) ' ' num2str(SmoothingIterations) ' ' OutDir '/HeadModel/m2m_' Subject '/SkinSmoothed.surf.gii']);
system(['wb_command -set-structure ' OutDir '/HeadModel/m2m_' Subject '/SkinSmoothed.surf.gii CORTEX_LEFT']);
system(['wb_command -set-structure ' OutDir '/HeadModel/m2m_' Subject '/Skin.surf.gii CORTEX_LEFT']);

% now, create a distance matrix; this is a vertex x vertex matrix
% summarizing the distance in geodesic space between all points on the
% head. Note: this can take awhile... 

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


