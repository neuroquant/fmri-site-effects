run('env.m'); 

dataset1 = load([GEN_DATAPATH filesep '3749595figshare/' filesep 'empirical_HCP']); 
dataset2 = load([GEN_DATAPATH filesep '3749595figshare/' filesep 'empirical_Paris']); 

% Create data matrix Y and covariate matrix X;
p = length(dataset1.TS) + length(dataset2.TS); 
n_rois = size(dataset1.TS{1},1);
n_volumes = min([size(dataset1.TS{1},2) size(dataset2.TS{1},2)]); 

Y = zeros(p,n_rois,n_volumes); 
X = zeros(p,1,n_volumes);
n_sites = 2;
standardize_rois = true;
for vol_no = 1:n_volumes;
	
	tmp_Y = zeros(p,n_rois);
		
	p1 = length(dataset1.TS);	
	for subj_no = 1:p1;
		tmp_Y(subj_no,:) = dataset1.TS{subj_no}(:,vol_no);
	end
	p2 = length(dataset2.TS);	
	for subj_no = 1:p2
		tmp_Y(p1+subj_no,:) = dataset2.TS{subj_no}(:,vol_no); 
	end
	
	tmp_X(1:p1,1) = 1; 
	tmp_X(p1+1:p,1) = 2;
	
	if(standardize_rois)
		std_y = std(tmp_Y); 
		Y(:,:,vol_no) = bsxfun(@rdivide,tmp_Y,std_y);
	else
		Y(:,:,vol_no) = tmp_Y;
	end
	X(:,1:n_sites,vol_no) = dummyvar(tmp_X);
	
	clear tmp_Y tmp_X;
end


% Correction per volume
results = {};
mu = zeros(n_volumes,n_rois); 
B = zeros(n_volumes,2,n_rois); 
Yres = zeros(size(Y)); 
U = zeros(size(Y)); 
for vol_no=1:n_volumes;
	results{vol_no} = projpca(Y(:,:,vol_no),X(:,1:2,vol_no)); 
	B(vol_no,:,:) = results{vol_no}.beta;
	mu(vol_no,:) = results{vol_no}.mu;	
	Yres(:,:,vol_no) = results{vol_no}.Y;
	U(:,:,vol_no) = results{vol_no}.U;
end

figure;
subplot(3,2,1); 
imagesc(corr(squeeze(Y(1,:,:))')); axis image equal; 
title('Original Site 1')
subplot(3,2,2); 
imagesc(corr(squeeze(Y(end,:,:))')); axis image equal; 
title('Original Site 2')
subplot(3,2,3); 
imagesc(corr(squeeze(Yres(1,:,:))')); axis image equal; 
title('Res + mean Site 1')
subplot(3,2,4); 
imagesc(corr(squeeze(Yres(end,:,:))')); axis image equal; 
title('Res + mean Site 2')
subplot(3,2,5); 
imagesc(corr(squeeze(U(1,:,:))')); axis image equal; 
title('Res Site 1')
subplot(3,2,6); 
imagesc(corr(squeeze(U(end,:,:))')); axis image equal; 
title('Res Site 2')

% Plot mean site effects
figure;
roi_idx = 2
subplot(3,2,1); 
plot(squeeze(Y(1,roi_idx,:))); 
title('Original Site 1')
subplot(3,2,2); 
plot(squeeze(Y(end,roi_idx,:))); 
title('Original Site 2')
subplot(3,2,3); 
plot(squeeze(Yres(1,roi_idx,:))); 
title('Res + mean Site 1')
subplot(3,2,4); 
plot(squeeze(Yres(end,roi_idx,:))); 
title('Res + mean Site 2')
subplot(3,2,5); 
plot(squeeze(U(1,roi_idx,:))); 
title('Res Site 1')
subplot(3,2,6); 
plot(squeeze(U(end,roi_idx,:))); 
title('Res Site 2')


% % Shared correction over all volumes
% % This approach is not recommended.
% % When rois are scaled to have unit variance, this also does nothing.
% results = projpca(reshape(Y,[p*n_volumes n_rois]),reshape(X,[p*n_volumes n_sites]));
% mu = results.mu;
% B = results.beta;
% Yres = reshape(results.Y,[p n_rois n_volumes]);
% U = reshape(results.U,[p n_rois n_volumes]);
%
% figure;
% subplot(3,2,1);
% imagesc(corr(squeeze(Y(1,:,:))')); axis image equal;
% title('Original Site 1')
% subplot(3,2,2);
% imagesc(corr(squeeze(Y(end,:,:))')); axis image equal;
% title('Original Site 2')
% subplot(3,2,3);
% imagesc(corr(squeeze(Yres(1,:,:))')); axis image equal;
% title('Res + mean Site 1')
% subplot(3,2,4);
% imagesc(corr(squeeze(Yres(end,:,:))')); axis image equal;
% title('Res + mean Site 2')
% subplot(3,2,5);
% imagesc(corr(squeeze(U(1,:,:))')); axis image equal;
% title('Res Site 1')
% subplot(3,2,6);
% imagesc(corr(squeeze(U(end,:,:))')); axis image equal;
% title('Res Site 2')
%
% % Plot mean site effects
% figure;
% roi_idx = 2
% subplot(3,2,1);
% plot(squeeze(Y(1,roi_idx,:)));
% title('Original Site 1')
% subplot(3,2,2);
% plot(squeeze(Y(end,roi_idx,:)));
% title('Original Site 2')
% subplot(3,2,3);
% plot(squeeze(Yres(1,roi_idx,:)));
% title('Res + mean Site 1')
% subplot(3,2,4);
% plot(squeeze(Yres(end,roi_idx,:)));
% title('Res + mean Site 2')
% subplot(3,2,5);
% plot(squeeze(U(1,roi_idx,:)));
% title('Res Site 1')
% subplot(3,2,6);
% plot(squeeze(U(end,roi_idx,:)));
% title('Res Site 2')
