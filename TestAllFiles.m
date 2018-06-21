classdef TestAllFiles < matlab.unittest.TestCase
    % To run:
    %  result = runtests('TestAllFiles')
    properties(TestParameter)
        file = getListOfFiles;
    end
    
    properties
        % These files require hardware, take too long to test,
        % are functions, or are templates
        toSkip = {'vittest.m','loopback.m','plutoFreqOffset.m',...
            'plutoLoopback.m','pluto_timing_offset.m',...
            'template_rt.m','transmitoffset.m','captureExample.m',...
            'template_1.m','template_2.m','template_3.m',...
            'cyclic.m','chanEQLMSDD.m','interpControl.m',...
            'interpFilter.m','loopFilter.m','transmitrepeat.m',...
            'transmitzeros.m','zcTED.m'}
    end
    
    methods(Test)
        function testFileDoesNotError(testCase, file)
            if ~any(strcmp(testCase.toSkip,file))
                fprintf('Testing: %s\n',file);
                run(fullfile('matlab', file));
                close all;
            end
        end
    end
end

function files = getListOfFiles
list    = dir(fullfile('matlab', '*.m'));
files = {list.name};
end
