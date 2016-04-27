opts.oversamp=2;
opts.N=16;
opts.k=2;
opts.num_disks=30;
opts.xsupp_msk=5;
%
clear intdim ecnt dindex j rand_seed Tv Nv x0 s0
%
mincos = .99;
randmin =1;
randmax =500;
%Build examples
ecnt = 0;
ftang=figure; set(ftang,'position',[100,1500,300,300]);

for j=randmin:randmax
    %PH_1 sets the random number generator through this variable
opts.seed=j;
[x0,s0]=PH_1(opts);
u=abs(fft0(ref_image));
algopts.support_mask=s0;
%Compute the tangent and normal spaces to the amplitude torus defined by x0
[Tv Nv] = torus_tansp(x0);
%Now we compare to the support constraint
fs0=s0(:);
szt=size(Tv);
fTv=reshape(Tv,[szt(1)*szt(2),szt(3)]);
[u d v]=svd(fTv((fs0 == 1),:));
dd=diag(d);
dindex=(dd>mincos);
%
ecnt=ecnt+1;
intdim(ecnt,1)=j;
intdim(ecnt,2)=sum(dindex);
        plot(dd(dindex))
        hold on
end
%Plot the histogram of intersection dimensions
hold off
figure, hist(intdim(:,2),20)
title(sprintf('k = %d, rmin = %d, rmax = %d, max-angle = %g, t-dim= %d',opts.k,randmin,...
    randmax, acos(mincos),length(x0(:))))
%Need to create a o/n basis for the support subspace to compare it to the
%tangent space to the torus.