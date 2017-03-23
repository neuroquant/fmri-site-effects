MATLABLIB = fullfile(getenv('PI_HOME'),'software','matlab','matlab-library');
if(~exist(MATLABLIB,'dir'))
	warning('Edit setup.m and set `MATLABLIB` to path to matlab-library');
	disp('To get matlab-library'); 
	disp('git clone git@github.com:TheEtkinLab/matlab-library.git'); 
	disp('Then follow instructions within matlab-library'); 
end
addpath(genpath(fullfile(MATLABLIB,'lib','ArgUtils'))); 
addpath(genpath(fullfile(MATLABLIB,'lib','BrewerMap'))); 
addpath(genpath(fullfile(MATLABLIB,'lib','gramm'))); 
addpath(genpath(fullfile(MATLABLIB,'test','matlab-xunit'))); 


if(exist('ggmClass') & exist('GGM','class'))
	warning('off')
	addpath('ggmClass');
	addpath('ggmClass/external'); 
	run('ggmClass/setup.m');
	warning('on');
else
	warning('ggmClass package is missing'); 
	disp('Attempting to download submodule ...'); 
	system('git submodule update --init ggmClass'); 
	warning('off')
	run('ggmClass/setup.m');
	warning('on');
end

if(exist('connectivity-diagnostics'))
	addpath('connectivity-diagnostics')
	if(exist('connectivity-diagnostics/startup.m'))
		run('connectivity-diagnostics/startup.m');
	end	
else
	warning('connectivity-diagnostics package is missing')
end