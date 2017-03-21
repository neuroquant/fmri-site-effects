MATLABLIB = fullfile(getenv('PI_HOME'),'software','matlab','matlab-library');
addpath(genpath(fullfile(MATLABLIB,'lib','ArgUtils'))); 
addpath(genpath(fullfile(MATLABLIB,'lib','BrewerMap'))); 
addpath(genpath(fullfile(MATLABLIB,'lib','gramm'))); 
addpath(genpath(fullfile(MATLABLIB,'test','matlab-xunit'))); 


if(exist('ggmClass'))
	warning('off')
	addpath('ggmClass');
	addpath('ggmClass/external'); 
	run('ggmClass/setup.m');
	warning('on');
else
	warning('ggmClass package is missing')
end

if(exist('connectivity-diagnostics'))
	addpath('connectivity-diagnostics')
	if(exist('connectivity-diagnostics/startup.m'))
		run('connectivity-diagnostics/startup.m');
	end	
else
	warning('connectivity-diagnostics package is missing')
end