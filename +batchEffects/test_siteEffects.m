%test_siteEffects.m
function test_suite = test_siteEffects
	buildFunctionHandleTestSuite(localfunctions);
	
end


function test_siteEffects_unit_diagonal
	
	p = 10; 
	n = 5; 
	
	Apool = eye(p); 
	Aind = repmat(Apool,[1 1 n]);
	
	% Check within_group_error
	metrics = siteEffects.within_group_error(...
						'Individual', Aind,...
						'Pooled', Apool ...
						);
						
	assert(any(metrics>=eps)==0, ...
			'siteEffects.within_group_error() fails null unit-diagonal test') 				
end