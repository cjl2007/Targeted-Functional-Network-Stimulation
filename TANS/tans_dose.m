function [OnTarget,Penalty,HotSpotSize] = tans_dose(magnE,VertexSurfaceArea,DiDt,AbsoluteThreshold,MinHotSpotSize,TargetNetwork,AvoidanceRegion,FunctionalNetworks,NetworkLabels,OutDir,Paths)
% cjl; cjl2007@med.cornell.edu;
%
% The inputs to the function are:
% magnE: A CIFTI file containing the e-field produced by the coil position prescribd by tans_optimize (can be either a string path or a structure array).
% VertexSurfaceArea: A CIFTI file containing the vertex surface areas (can be either a string path or a structure array). 
% DiDt: A range of stimulation intensities (the speed of variation of the current throughout the stimulating coil, in units of A/us, numeric). 
% AbsoluteThreshold: This value (in V/m units) represents an assumed neural activation threshold (numeric). 
% MinHotSpotSize: The min. size (in surface area, mm2) of the E-field hotspot. 
% TargetNetwork: A CIFTI file containing the functional network of interest (structure array). Non-zero values in TargetNetwork.data are considered target network vertices.
% AvoidanceRegion: A CIFTI file indexing the brain regions or networks of no interest (structure array). Non-zero values in AvoidanceRegion.data are considered non-target network vertices. 
% FunctionalNetworks: A CIFTI file indexing all functional networks.
% NetworkLabels: A .xls database describing the Names and Colors of each functional networks.
% OutDir: Path to the output folder (string).
% Paths: Paths to folders that must be added to Matlab search path (cell array of strings). 

% add some directories 
% to the search path;
for i = 1:length(Paths)
addpath(genpath(Paths{i})); % 
end

% read in the standardized network file;
NetworkLabels = dataset('XLS',NetworkLabels);

% read in 
% the E-field;
if ischar(magnE)
    magnE = ft_read_cifti_mod(magnE);
end

% read in the
% surface area;
if ischar(VertexSurfaceArea)
    VertexSurfaceArea = ft_read_cifti_mod(VertexSurfaceArea);
end

% sweep intensity levels
for i = 1:length(DiDt)
    magnE.data(:,i) = (magnE.data(:,1) / (1*1e6)) * DiDt(i); 
end

% write out the cifti file; each column == a different stimulation intensity level
ft_write_cifti_mod([OutDir '/magnE_BestCoilCenter+BestOrientation_AllStimulationItensities.dtseries.nii'],magnE);
system(['echo ' num2str(DiDt/1e6) ' > ' OutDir '/DiDt.txt']); % write out the intensities used;

% if no avoidance
% region is specified;
if isempty(AvoidanceRegion)
    AvoidanceRegion = TargetNetwork; % preallocate
    AvoidanceRegion.data = zeros(size(AvoidanceRegion.data,1));
end

% preallocate;
OnTarget = zeros(size(magnE.data,2),1); % "OnTaget" variable (% of E-field hotspot that contains target network vertices);
Penalty = zeros(size(magnE.data,2),1); % "Penalty" variable (% of E-field hotspot that contains avoidance region / network vertices);
HotSpotSize = zeros(size(magnE.data,2),1);

% preallocate;
O = TargetNetwork;
O.data = zeros(size(O.data,1),length(DiDt)); % clear contents of the file;

% log the functional networks inside the blank cifti multiple times;
O.data(1:59412,1:length(DiDt)) = repmat(FunctionalNetworks.data(1:59412,1),[1 length(DiDt)]);
 
% sweep the range of
% stimulation intensities
for t = 1:size(magnE.data,2)
    
    HotSpot = magnE.data(1:59412,t) >= AbsoluteThreshold; % this is the hotspot (suprathreshold portion of e-field);
    HotSpotSize(t) = sum(VertexSurfaceArea.data(HotSpot)); % this is the total surface area of the hotspot 
    
    % calculate proportion of the suprathreshold E-field that is on-target
    OnTarget(t,1) = sum(VertexSurfaceArea.data(HotSpot&TargetNetwork.data(1:59412,1)==1)) / sum(VertexSurfaceArea.data(HotSpot));
    Penalty(t,1) = sum(VertexSurfaceArea.data(HotSpot&AvoidanceRegion.data(1:59412,1)==1)) / sum(VertexSurfaceArea.data(HotSpot));

    % discard functional networks
    % outside of E-field hotspot;
    O.data(HotSpot==0,t)=0;
    
end

% write out the stimulated networks;
ft_write_cifti_mod([OutDir '/StimulatedNetworks'],O);

% write out the first community.
system(['echo ' char(NetworkLabels.Network(1)) ' > ' OutDir '/LabelListFile.txt']);
system(['echo 1 ' num2str(NetworkLabels.R(1)) ' ' num2str(NetworkLabels.G(1)) ' ' num2str(NetworkLabels.B(1)) ' 255 >> ' OutDir '/LabelListFile.txt ']);

% sweep through 
% the other networks;
for i = 2:size(NetworkLabels,1)
    system(['echo ' char(NetworkLabels.Network(i)) ' >> ' OutDir '/LabelListFile.txt']);
    system(['echo ' num2str(i) ' ' num2str(NetworkLabels.R(i)) ' ' num2str(NetworkLabels.G(i)) ' ' num2str(NetworkLabels.B(i)) ' 255 >> ' OutDir '/LabelListFile.txt ']);  
end

% make dlabel cifti file summarizing which networks were stimulated;
system(['wb_command -cifti-label-import ' OutDir '/StimulatedNetworks.dtseries.nii ' OutDir '/LabelListFile.txt ' OutDir '/StimulatedNetworks_AbsThreshold_' num2str(AbsoluteThreshold) 'Vm_MinHotSpotSize_' num2str(MinHotSpotSize) 'mm2.dlabel.nii -discard-others']);
system(['rm ' OutDir '/StimulatedNetworks*.dtseries.nii']); % clear some intermediate files;
system(['rm ' OutDir '/LabelListFile.txt']); % clear some intermediate files;

PenalizedOnTarget = OnTarget - Penalty;
PenalizedOnTarget(HotSpotSize < MinHotSpotSize) = 0; 

% find the best dose;
Idx = find(PenalizedOnTarget == max(PenalizedOnTarget));
magnE.data = magnE.data(:,Idx);

% write out the cifti file; 
ft_write_cifti_mod([OutDir '/magnE_BestCoilCenter+BestOrientation+BestDose.dtseries.nii'],magnE);
system(['echo ' num2str(DiDt(Idx)/1e6) ' > ' OutDir '/BestDose.txt']); % write out the best dose.

% make some figures;

H = figure; % prellocate parent figure
set(H,'position',[1 1 350 225]); hold;

% plot the results;
plot(PenalizedOnTarget * 100,'Color','k','LineWidth',1); hold;

% make it "pretty";
xlim([0 length(DiDt)]);
ylim([0 100]);
yticks(0:20:100);
xlabel('dI/dt (A/\mus)');
ylabel('% On-Target')
Idx = 1:round(length(DiDt)/20):length(DiDt);
xticks(Idx);
xticklabels(round(DiDt(Idx)/1e6));
xtickangle(90);
set(gca,'FontName','Arial','FontSize',12,'TickLength',[0 0]);
saveas(gcf,[OutDir '/OnTarget_AbsThreshold_' num2str(AbsoluteThreshold) 'Vm_MinHotSpotSize_' num2str(MinHotSpotSize) 'mm2.pdf']); 
close all;

H = figure; % prellocate parent figure
set(H,'position',[1 1 350 225]); hold;

% plot the results;
OnTarget(isnan(OnTarget))=0;
scatter(HotSpotSize,OnTarget * 100,[],DiDt/1e6,'filled'); 
colormap(jet);
h = colorbar; 
h.Label.String = 'dI/dt (A/\mus)';

% make it "pretty";
ylim([0 100]);
yticks(0:20:100);
xlabel('Hotspot size (mm2)');
ylabel('% On-Target')
set(gca,'FontName','Arial','FontSize',12,'TickLength',[0 0]);
saveas(gcf,[OutDir '/OnTarget_vs_HotspotSize_Curve_AbsThreshold_' num2str(AbsoluteThreshold) 'Vm.pdf']); 
close all;

end
