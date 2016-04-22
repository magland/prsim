%Options to define the phantom
%clear all
opts.oversamp=2;
opts.N=16;
opts.k=2;
opts.num_disks=30;
%Give the variance for the initial choice of phases. 
%For completely random choice set to opts.N
opts.sigma=4;%opts.N;
%
%Options to control the auxiliary conditions used
algopts.positivity=0;
algopts.support=1;
%To regularize the support condition
algopts.eps_reg=0;
mincos =.98;
rand_seed = 762
%How many iterations should we try?
num_iterations=300;
%Build examples
%Choose an algorithm
ALG=@ALG_NEWT;
Tvflag =0;
%
%ALG=@ALG_BA_HIO;
%ALG=@ALG_AB_HIO;
%
%If two_step = 1, then we run a second algorithm starting with the output
%of the first one.
two_step = 0;
%
algopts.alg=ALG;
if strcmp(func2str(ALG),'ALG_NEWT')
    Tvflag =1;
end
%Set the random number generator
%
rng(rand_seed);
[ref_image,support_mask]=PH_1(opts);
algopts.support_mask=support_mask;
u=abs(fftn(ref_image));
init_image = 0;
%If we have a good intial image to try
%init_image = recon;
%
%Randomize phases to start test
%Need to use fft0 to get the initial guess to be a real image
%We randomize around the known phases of the initial image
if (init_image == 0)
    u0=fft0(ref_image);
    sigma=opts.sigma/opts.N %How random is the initial data sigma = 1 is completely random
                   %We divide by opts.N to control the L^2 error
    dth=2*pi*sigma*(rand(size(u0))-0.5*ones(size(u0)));

    init_image_hat=u0.*exp(dth*i);
    init_image=ifft0(init_image_hat);
    init_err=norm(ref_image(:)-init_image(:))/norm(ref_image(:));
    %
    figure, hist(2*dth(:)/pi,20)
    title(sprintf('Initial phase errors mod pi/2, err = %g',init_err))
end
%
%
DD=struct;
DD.x=init_image;
if (exist('fA','var')) clear fA; end;
fA=figure;
alg_resids=[];
errors=[];
ecnt =0;
 for it=1:num_iterations
            DD.itnum=it;
            [DD,recon,info]=ALG(DD,u,algopts);
            %The flattened tangent space
            if (Tvflag == 1)
                fs0=support_mask(:);
                szt=size(info.Tv);
                fTv=reshape(info.Tv,[szt(1)*szt(2),szt(3)]);
                %We compute the intersection between the linear subspace
                %underlying the (affine) tangent space to the support and
                %the subspace defined by the support condition.
                [uu d vv]=svd(fTv((fs0 == 1),:));
                dd=diag(d);
                dindex=(dd>mincos);
                %
                ecnt=ecnt+1;
                intdim(ecnt,1)=it;
                intdim(ecnt,2)=sum(dindex);
            end
            alg_resids(it)=info.resid;
            tA=tic;
            
            %if (toc(tA)>1)||(it==1)||(it==num_iterations)
                recon=register_to_reference(recon,ref_image);
                errors(it)=compute_residual(recon,ref_image);
                update_single_run_figure(fA,ref_image,init_image,recon,alg_resids,...
                    errors,algopts,opts);
                %drawnow;
                %tA=tic;
%             else
%                 errors(it)=nan;
%             end;
end;
%We can use the ouput of the previous run as the input to a different algorithm
%
if two_step == 1
    %Options to control the auxiliary conditions used
    algopts2.positivity=1;
    algopts2.support=0;
    %To regularize the support condition
    algopts2.eps_reg=0%10^(-8);
    algopts2.support_mask=support_mask;
    %How many iterations should we try?
    num_iterations2=4;
    %Build examples
    %Choose an algorithm
    %ALG=@ALG_NEWT;
    Tvflag2 =0;
    %
    ALG2=@ALG_NEWT;
    algopts2.alg=ALG2;
     if strcmp(func2str(ALG2),'ALG_NEWT')
        Tvflag2 =1;
     end   
    %We initialize with the output of the previous algorithm.
    DD2=struct;
    DD2.x=recon;
    if (exist('fA2','var')) clear fA2; end;
    fA2=figure;
    alg2_resids=[];
    errors2=[];
    ecnt2 =0;
     for it=1:num_iterations2
                DD2.itnum=it;
                [DD2,recon2,info2]=ALG2(DD2,u,algopts2);
                %The flattened tangent space
                if (Tvflag2 == 1)
                    fs0=support_mask(:);
                    szt=size(info2.Tv);
                    fTv=reshape(info2.Tv,[szt(1)*szt(2),szt(3)]);
                    [uu d vv]=svd(fTv((fs0 == 1),:));
                    dd=diag(d);
                    dindex=(dd>mincos);
                    %
                    ecnt2=ecnt2+1;
                    intdim2(ecnt2,1)=it;
                    intdim(ecnt2,2)=sum(dindex);
                end
                alg2_resids(it)=info2.resid;
                tA=tic;

                %if (toc(tA)>1)||(it==1)||(it==num_iterations)
                    recon2=register_to_reference(recon2,ref_image);
                    errors2(it)=compute_residual(recon2,ref_image);
                    update_single_run_figure(fA2,ref_image,recon,recon2,alg2_resids,...
                        errors2,algopts2,opts);
                    %drawnow;
                    %tA=tic;
    %             else
    %                 errors(it)=nan;
    %             end;
    end;
end