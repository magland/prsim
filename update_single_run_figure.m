function update_single_run_figure(fA,ref_image,init_image,recon,alg_resids,errors,algopts,opts)
figure(fA);

resid=abs(ref_image-recon);
subplot(2,3,1:3);
show_image(cat(1,ref_image,recon,resid));
title(sprintf(' N= %d, sigma = %d, pos = %d, supp = %d,  eps_reg = %g,   ALG= %s',...
    opts.N,opts.sigma,algopts.positivity,algopts.support,algopts.eps_reg,func2str(algopts.alg)));

subplot(2,3,4);
show_image(init_image)
title('Init Image')
subplot(2,3,5);
semilogy(1:length(alg_resids),alg_resids);
title('Residuals')
inds0=find(~isnan(errors));
subplot(2,3,6);
semilogy(inds0,errors(inds0));
title('Errors')
end
