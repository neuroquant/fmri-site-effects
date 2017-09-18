classdef batchEffects <  matlab.mixin.Copyable
%BATCHEFFECTS - General functions to diagnose and quantify session, site or other batch effects for matrix data. 
% 
%FUNCTIONS: 
% 	
% 	- siteEffects.test_retest_mccc
% 	  metrics = siteEffects.test_retest_mccc(X)
%
% 	- siteEffects.within_group_error 	
% 	  metrics = siteEffects.within_group_error('Individual', Aind, 'Pooled', Apooled); 
% 
% 
	methods
		
		function disp(self)
			
			fprintf(help('batchEffects')); 
			
		end

	end

	methods(Static)

        function metrics = test_retest_mccc(varargin)
        % TEST_RETEST_MCCC - Calculates the multivariate concordance correlation coefficient to assess reliability of $p$ features between two raters using n i.i.d samples (for e.g. subjects) 
        %  
        % USAGE: 
        %   metrics = batchEffects.test_retest_mccc(X)
        % 
        % INPUTS: 
        %   X - is a n_samples x n_features x 2 raters data matrix
        % 
        % 
        
            
        end
        

		function metrics = within_group_variation(varargin)
		%WITHIN_GROUP_VARIATION - Calculates matrix norm between individual matrices and the pooled matrix.
		% 
		% USAGE: 
		% 	metrics = batchEffects.within_group_variation('Individual', Aind, 'Pooled', Apooled); 
		% 
		% INPUTS: 'Param1', Value1, 'Param2', Value2, ... 	
		% 	'Individual', 	Aind is a p x p x n matrix such as n subject level correlation matrices of size p x p
		% 	'Pooled', 		Apooled is a single p x p matrix such as group correlation matrix
		% 	options 		Not yet supported. 
		% 	options.norm 	'fro','Inf', 1, 2, 'geo'
	
			narginchk(4,4);
			warning('Input arguments are not being checked. TBA'); 

			Aind 	= varargin{2}; 
			Apooled = varargin{4}; 
	
			if(ndims(Aind)==2)
				n=1;
			elseif(ndims(Aind)==3)
				n = size(Aind,3);
			else
				error('Expecting p x p x n matrix for `Individual` parameter ');
			end

			metrics = zeros(1,n); 
			for matrix_no = 1:n
				metrics(matrix_no) = norm(Aind(:,:,matrix_no)-Apooled,'fro'); 		
			end
		end	

	end
end