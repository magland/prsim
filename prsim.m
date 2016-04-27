function results=prsim(ALG,PHANTOM,algopts,phopts,simopts)

if nargin<1, test_prsim; return; end;

if (isfield(simopts,'seed')) rng(simopts.seed); end;

num_trials=simopts.num_trials;
num_restarts=simopts.num_restarts;
num_iterations=simopts.num_iterations;

results.errors=zeros(num_restarts,num_trials);
for tt=1:num_trials
    [ref_image,support_mask]=PHANTOM(phopts);
    algopts.support_mask=support_mask;
    u=abs(fft0(ref_image));
    for rr=1:num_restarts
        init_image_hat=u.*exp(2*pi*i*rand(size(u)));
        init_image=ifft0(init_image_hat);
        DD=struct;
        DD.x=init_image;
        if (exist('fA','var')) close(fA); end;
        fA=figure;
        tA=tic;
        alg_resids=[];
        errors=[];
        for it=1:num_iterations
            DD.itnum=it;
            [DD,recon,info]=ALG(DD,u,algopts);
            alg_resids(it)=info.resid;
            if (toc(tA)>1)||(it==1)||(it==num_iterations)
                recon=register_to_reference(recon,ref_image);
                errors(it)=compute_residual(recon,ref_image);
                update_single_run_figure(fA,ref_image,init_image,recon,alg_resids,errors,algopts,phopts);
                title(sprintf('trial %d, restart %d',tt,rr));
                drawnow;
                tA=tic;
            else
                errors(it)=nan;
            end;
        end;
        pause(0.2);
        results.errors(rr,tt)=errors(end);
        results.alg_resids{rr,tt}=alg_resids;
    end;
end;

figure; hist(results.errors,100)
title(sprintf('Errors, Algorithm %s, positivity %d, support %d',...
info.name,algopts.positivity,algopts.support));

figure;
for tt=1:num_trials
for rr=1:num_restarts
    tmp=results.alg_resids{rr,tt};
    semilogy(1:length(tmp),tmp); hold on;
end;
title(sprintf('Algorithm %s, positivity %d, support %d',info.name,algopts.positivity,algopts.support));
end;

end

function update_single_run_figure(fA,ref_image,recon,alg_resids,errors)
figure(fA);

resid=abs(ref_image-recon);
subplot(2,3,1:3);
show_image(cat(1,ref_image,recon,resid));

subplot(2,3,4);
semilogy(1:length(alg_resids),alg_resids);

inds0=find(~isnan(errors));
subplot(2,3,5);
semilogy(inds0,errors(inds0));
end

function show_image(X)
imagesc(X'); colormap('gray');
end

