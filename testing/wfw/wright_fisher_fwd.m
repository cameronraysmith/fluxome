function [Zs, Xs, Pis, log_Ps] = wright_fisher_fwd(N,T,nSim,theta_f,theta_h,theta_z0,...
    theta_g,bin_expr_flag,fwd_sigma,seed,verbose)

Zs = cell(1,nSim);
Xs = cell(1,nSim);
Pis = cell(1,nSim);
log_Ps = zeros(1,nSim);
Dx = length(theta_f);
Dz = size(theta_g,1);

rng(seed);

for cSim = 1:nSim

    if verbose
        cSim
    end
    Z = zeros(T,N,Dz);
    X = zeros(T,N,Dx);
    Pi = zeros(T,N);
    for i = 1:Dz
        Z(1,:,i) = (rand(1,N) < theta_z0(i));
        log_Ps(cSim) = log_Ps(cSim) + sum(Z(1,:,i)*log(theta_z0(i)) + (1-Z(1,:,i))*log(1-theta_z0(i)));
    end

    preds = squeeze(Z(1,:,:)) * theta_g;
    if bin_expr_flag==1
        preds = 1 ./ (1 + exp(-preds));
        for i = 1:Dx
            X(1,:,i) = (rand(1,N) < preds(:,i)');
            log_Ps(cSim) = log_Ps(cSim) + sum(X(1,:,i).*log(preds(:,i)') + (1-X(1,:,i)).*log(1-preds(:,i)'));
        end
    else
        mat = mvnrnd(preds,fwd_sigma*eye(Dx));
        X(1,:,:) = mat;
        [dum vec] = mvnpdf_log(squeeze(X(1,:,:)),preds,fwd_sigma*eye(Dx));
        log_Ps(cSim) = log_Ps(cSim) + sum(vec);
    end

    for t = 2:T
        fs = zeros(1,N);
        for i = 1:Dx
            fs = fs + X(t-1,:,i) * theta_f(i);
        end
        fs = exp(fs);
        fs = fs ./ sum(fs);
        fs_ = cumsum(fs);
        for n = 1:N
            idx = find(rand<=fs_);
            idx = idx(1);
            Pi(t,n) = idx;            
            log_Ps(cSim) = log_Ps(cSim) + log(fs(idx));
            for i = 1:Dz
                if (rand >= theta_h(i))
                    Z(t,n,i) = Z(t-1,idx,i);
                    log_Ps(cSim) = log_Ps(cSim) + log(1-theta_h(i));
                else
                    Z(t,n,i) = 1 - Z(t-1,idx,i);
                    log_Ps(cSim) = log_Ps(cSim) + log(theta_h(i));
                end
            end
        end
        preds = squeeze(Z(t,:,:)) * theta_g;
        if bin_expr_flag==1
            preds = 1 ./ (1 + exp(-preds));
            for i = 1:Dx
                X(t,:,i) = (rand(1,N) < preds(:,i)');
                log_Ps(cSim) = log_Ps(cSim) + sum(X(t,:,i).*log(preds(:,i)') + (1-X(t,:,i)).*log(1-preds(:,i)'));
            end
        else
            mat = mvnrnd(preds,fwd_sigma*eye(Dx));
            X(t,:,:) = mat;
            [dum vec] = mvnpdf_log(squeeze(X(t,:,:)),preds,fwd_sigma*eye(Dx));
            log_Ps(cSim) = log_Ps(cSim) + sum(vec);
        end
    end
    Zs{cSim} = Z;
    Xs{cSim} = X;
    Pis{cSim} = Pi;    
end
