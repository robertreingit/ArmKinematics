%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get_marker
% Simple extractor function to obtain the 6D marker from the from the 
% data matrix R^(m x no_markers*6), 
%       m = number of frames
%       no_markers = number of markers
% INPUT:
% data  = marker data
% idx   = number of marker (>=1)
% OUTPUT:
% m = marker data e R^mx6 containing the 
%       position data => columns 1-3
%       attitude data => columns 4-6
% SIDEEFFECTS:
% None.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = get_marker(data,idx)

    m = data(:,(idx-1)*6 + (1:6));

end
