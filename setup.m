MATLABLIB = fullfile(getenv('PI_HOME'),'software','matlab','matlab-library');
addpath(genpath(MATLABLIB)); 

if(exist('ggmClass'))
	addpath(ggmClass);
	addpath(ggmClass/external); 
	run('ggmClass/setup.m');
end