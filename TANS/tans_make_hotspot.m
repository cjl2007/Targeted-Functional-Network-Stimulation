function [HotSpot] = tans_make_hotspot(E,Thresholds)

% preallocate the E-field hot spot variable;
HotSpot = zeros(size(E.data,1),length(Thresholds));

% sweep the E-field thresholds
for i = 1:length(Thresholds)
    
    % this is the E-field "hotspot" at threshold "i";
    HotSpot(E.data > prctile(E.data,Thresholds(i)),i) = true;
    
end

% concatenate across thresholds;
HotSpot = [sum(HotSpot,2) HotSpot];

end