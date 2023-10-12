function train_wright_fisher_model_wasserstein_gt

%%% evaluates wasserstein cost for ground truth objective 
%%% NB doesn't train a new model

rng(100);
N = 25; % size of population
Dz = 4; % dimensionality of Z (# variants)
Dx = 2; % dimensionality of X (# genes) 
T = 10;  % # time-points
nSim = 20; % # fwd simulations
theta_f = [1 -1]';  % log relative fitness of gene
theta_h = 0.05 * ones(Dz,1); % mutation rate
theta_z0 = 0.5 * ones(Dz,1); % initial probability of variants
theta_g_mask = [repmat([1 0],[2 1]) ; repmat([0 1],[2 1])]; % mask for g-p map
theta_g = [1 0 ; -1 0 ; 0 1 ; 0 -1]; % true g-p map
bin_expr_flag = 0; % binary expression flag
fwd_sigma = 0.01; % gaussian sig value for fwd model
verbose = 0; % verbosity
    
nSamp = 5;
max_mu = 3;
max_sig = 2;
nEpoch = 1;

sig = 0;
theta_f_mu = theta_f';
theta_g_mu = theta_g;
alpha = 1; % weighting of transcriptome loss
emd_norm = 2;
pFlip = 0.01;  % for variant matching
inv_temp = 1; % temperature

sd = 10;
[Zs, Xs, Pis, log_Ps] = wright_fisher_fwd(N,T,nSim,theta_f,theta_h,theta_z0,theta_g,bin_expr_flag,fwd_sigma,sd,verbose);

theta_f_mus_all = cell(1,nEpoch);
theta_g_mus_all = cell(1,nEpoch);
sigs_all = zeros(1,nEpoch);
scs_all = cell(1,nEpoch);
Fs_all = zeros(nSamp,nEpoch);
emds_all = zeros(nSamp,nEpoch);
theta_f_samps_all = cell(1,nEpoch);
theta_g_samps_all = cell(1,nEpoch);
sims_Zs_all = cell(1,nEpoch);
sims_Xs_all = cell(1,nEpoch);

rng(100);

for cEpoch = 1:nEpoch
    
%     cEpoch
    
    theta_f_samps = repmat(theta_f_mu,[nSamp 1]);
    theta_f_samps = theta_f_samps + (sig^2) * randn(size(theta_f_samps));
    theta_g_samps = repmat(theta_g_mu,[1 1 nSamp]);
    theta_g_samps = theta_g_samps + (sig^2) * randn(size(theta_g_samps));
    for i = 1:nSamp
        theta_g_samps(:,:,i) = theta_g_samps(:,:,i) .* theta_g_mask;
    end
    
    sims_Zs = cell(nSamp,nSim);
    sims_Xs = cell(nSamp,nSim);
    
    for i = 1:nSamp        
        sd1 = (2^sd) * (3^cEpoch) * (5^i);
        [Zs1, Xs1, Pis1, log_Ps1] = wright_fisher_fwd(N,T,nSim,theta_f_samps(i,:)',theta_h,theta_z0,theta_g_samps(:,:,i),bin_expr_flag,fwd_sigma,sd1,verbose);
        for j = 1:nSim
            sims_Zs{i,j} = Zs1{j};
            sims_Xs{i,j} = Xs1{j};
        end
    end

    emds = zeros(nSamp,nSim,nSim);
    for i = 1:nSamp   
        [cEpoch i]
        for j = 1:nSim
            for k = 1:nSim
%                 [i j k]
                Zs1 = squeeze(sims_Zs{i,j}(end,:,:));
                Xs1 = squeeze(sims_Xs{i,j}(end,:,:));
                Zs2 = squeeze(Zs{k}(end,:,:));
                Xs2 = squeeze(Xs{k}(end,:,:));  
                
                [dum fval] = emd([Zs1*alpha Xs1], [Zs2*alpha Xs2], ones(N,1)./N, ones(N,1)./N, (@(V1,V2) norm(V1 - V2, emd_norm)));
                
%                 mean(Xs2)
                
                emds(i,j,k) = fval;

            end
        end
    end
    F_scs = zeros(nSamp,1);
    for n = 1:nSamp
        F_scs(n) = -emd2(squeeze(emds(n,:,:)));
    end
    emds_all(:,cEpoch) = -F_scs;
    F_scs1 = F_scs - max(F_scs);
    F_scs = exp(inv_temp * F_scs); 
    F_scs1 = exp(inv_temp * F_scs1); 
    
    theta_f_mus_all{cEpoch} = theta_f_mu;
    theta_g_mus_all{cEpoch} = theta_g_mu;
    sigs_all(cEpoch) = sig;
    scs_all{cEpoch} = emds;
    Fs_all(:,cEpoch) = F_scs;   
    theta_f_samps_all{cEpoch} = theta_f_samps;
    theta_g_samps_all{cEpoch} = theta_g_samps;
    sims_Zs_all{cEpoch} = sims_Zs;
    sims_Xs_all{cEpoch} = sims_Xs;
    
    sigs_all
    [mean(emds_all) ;  max(emds_all)]
    
    % update
    
%     if mod(cEpoch,2)==1
        ws = F_scs1 ./ sum(F_scs1);
        theta_f_mu = sum(theta_f_samps .* repmat(ws,[1 length(theta_f)]));
        ws_block = ones(size(theta_g_samps));
        for n = 1:nSamp
            ws_block(:,:,n) = ws(n);
        end
        theta_g_mu = sum(theta_g_samps .* ws_block,3);
%     else
        ws = F_scs1 ./ sum(F_scs1);
        vars = zeros(nSamp,1);
        for n = 1:nSamp
            mat = theta_g_samps(:,:,n);
            vars(n) = sum((theta_f_samps(n,:)-theta_f_mu).^2) + ...
                sum((mat(:)-theta_g_mu(:)).^2);
            vars(n) = vars(n) / (sum(theta_g_mask(:)) + length(theta_f));
        end
        sig = sqrt(sum(ws.*vars));
%     end

%     max(max(abs(theta_f_mu(:))),max(abs(theta_g_mu(:))))   
    
end

save('training_outputs_gt','theta_f_mus_all','theta_g_mus_all','sigs_all','scs_all',...
    'Fs_all','theta_f_samps_all','theta_g_samps_all',...
    'sims_Zs_all','sims_Xs_all','Zs', 'Xs', 'Pis', 'log_Ps', 'emds_all',...
    'theta_f_mu','theta_g_mu','sig','emds');
