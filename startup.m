run('setup.m');
[sysinfo sysval] = system('./update.sh'); 

if(sysinfo)
	warning('./update.sh failed');
	disp(sysval);
end