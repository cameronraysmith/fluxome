function plot_training_stats

load('training_outputs','theta_f_mus_all','theta_g_mus_all','sigs_all','scs_all',...
    'emds_all','theta_f_samps_all','theta_g_samps_all',...
    'sims_Zs_all','sims_Xs_all','Zs', 'Xs', 'Pis', 'log_Ps');

close all

figure(1);
plot(sigs_all,'g-','linewidth',2);

figure(2);
plot(mean(emds_all),'b-','linewidth',2); hold on;
plot(max(emds_all),'r-','linewidth',2);

best_f = [];
best_g = [];
rng(100);
Dz = 4; % dimensionality of Z (# variants)
Dx = 2; % dimensionality of X (# genes) 
theta_f = [1 -1]';  % log relative fitness of gene
theta_g_mask = [repmat([1 0],[2 1]) ; repmat([0 1],[2 1])]; % mask for g-p map
theta_g = [1 0 ; -1 0 ; 0 1 ; 0 -1]; % true g-p map

for i = 1:length(sigs_all)
    idx = find(emds_all(:,i)==min(emds_all(:,i)),1);
    best_f = [best_f ; theta_f_samps_all{i}(idx,:)];
    best_g = cat(3,best_g, theta_g_samps_all{i}(:,:,idx));
end

%%%
% plot

theta_f_ests = best_f';
theta_g_ests = best_g;

Dx = length(theta_f);
Dz = size(theta_g,1);
nEpoch = size(theta_f_ests,2)-1;
figure(3);
cols = {'k' 'r' 'b' 'g' 'c' 'm'};
for i = 1:Dx
    plot(0:nEpoch,theta_f_ests(i,:),[cols{i} '-'],'linewidth',1.5); hold on;
    plot([0 nEpoch],[theta_f(i) theta_f(i)],[cols{i} '--'],'linewidth',2); hold on;
end
% ylim([min(theta_f)-0.2, max(theta_f)+0.1]);

figure(4)
for ii = 1:Dx
    subplot(1,Dx,ii)
    idxs = find(theta_g_mask(:,ii)==1);
    for i = 1:length(idxs)
        vec = theta_g_ests(idxs(i),ii,:);
        plot(0:nEpoch,squeeze(vec),[cols{i} '-'],'linewidth',1.5); hold on;
        plot([0 nEpoch],[theta_g(idxs(i),ii) theta_g(idxs(i),ii)],[cols{i} '--'],'linewidth',2); hold on;
    end
end
% ylim([min(theta_g(:))-0.1, max(theta_g(:))+0.1]);

