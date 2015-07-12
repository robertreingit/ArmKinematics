classdef getMarkerTest < matlab.unittest.TestCase
   
    methods(Test)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % extractMarker
        % Get the data from marker at second position
        function extractMarker(testCase)
            
            data = repmat(1:12,[10,1]);
            testCase.verifyTrue(all(all(get_marker(data,2) == data(:,7:12))));
            
        end
        
    end
    
end