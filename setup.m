% Use matlab package manager to install requirements.txt
USE_MPM = false;

% Specify path to matlab-library
MATLABLIB = fullfile(getenv('PI_HOME'),'software','matlab','matlab-library');
USE_MATLAB_LIBRARY = true;

if(USE_MPM)
	if(exist('mpm'))
		disp('Invoking Matlab Package Manager')
	else
		disp('Installing Matlab Package Manager'); 
		disp('git clone https://github.com/mobeets/mpm.git external/mpm'); 
	end
	
	disp('Installing package dependencies to external/.')
	PYTHON_EXE = '/usr/bin/python';
	MPM_INSTALL_DIR = [pwd filesep 'external'];
	try
		mpm -r requirements.txt
		mpm_paths(MPM_INSTALL_DIR);		
	catch me 
		disp(me)
		warning('Could not install packages in requirements.txt. Please install them manually to external/ or add pre-existing packages to your path'); 
	end
	MATLABLIB = [pwd filesep 'external'];
end

if(USE_MATLAB_LIBRARY)	
	if(~exist(MATLABLIB,'dir'))
		warning('Edit setup.m and set `MATLABLIB` to path to matlab-library');
		% disp('To get matlab-library');
		% disp('git clone git@github.com:TheEtkinLab/matlab-library.git');
		% disp('Then follow instructions within matlab-library');	
	else
		addpath(genpath(fullfile(MATLABLIB,'lib','ArgUtils'))); 
		addpath(genpath(fullfile(MATLABLIB,'lib','BrewerMap'))); 
		addpath(genpath(fullfile(MATLABLIB,'lib','gramm'))); 
		addpath(genpath(fullfile(MATLABLIB,'test','matlab-xunit'))); 
	end
end

% Install submodule if missing
if(exist('ggmClass') & exist('ggmClass/setup.m','file'))
	warning('off')
	addpath('ggmClass');
	addpath('ggmClass/external'); 
	run('ggmClass/setup.m');
	warning('on');
else
	warning('ggmClass package is missing'); 
	disp('Attempting to download submodule ...'); 
	system('git submodule update --init ggmClass'); 
	warning('off');
	run('ggmClass/setup.m');
	warning('on');
end

% Install submodule if missing
if(exist('connectivity-diagnostics') & exist('connectivity-diagnostics/startup.m','file'))
	addpath('connectivity-diagnostics')
	if(exist('connectivity-diagnostics/startup.m'))
		run('connectivity-diagnostics/startup.m');
	end	
else
	warning('connectivity-diagnostics package is missing')
	disp('Attempting to download submodule ...'); 
	system('git submodule update --init connectivity-diagnostics'); 
	warning('off'); 
	run('connectivity-diagnostics/startup.m');
	warning('on');
end
