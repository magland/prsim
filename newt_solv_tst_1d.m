%This script is useful for studying the tangent space of a torus
%defined by the magnitude data of an FFT of a real sequence.
%
%If we set mnmx=0, then the routine will minimize the sum of the entries over
%indices in the complement of the "known" support of
%points in the affine tangent space to the amplitude torus intersected with
%the non-negative orthant. This is a linearized version of the phase retrieval
%problem: find a vector whose Fourier transform has given magnitudes with a 
%known support.
%
%If mnmx=1 then it maximizes this sum, thereby searching for points in the
%intersection of this tangent space with the interior of the non-negative
%orthant.
%
%
%Half the number of real parameters:
N=128;
%The scaled size of the perturbation in the phases:
sigma_0=.1;
%Random seed for creating the example
rand_seed=25;
%For computing intersection dimension
mincos =.99;
%
plt.singvec = 1;% Show plots of the singular vectors in the intersection of 
               % the tangent space.
%The rand_seed=6 produces an example where the tangent space to the torus
%is essentially transversal to the tangent space of the support subspace.
%The rand_seed=7 produces an example with a high order of tangency between
%these subspaces. For these examples
        % optsexam.k=4;
        % optsexam.num_disks=4;
        % rand_seed=6;
%zero padding within the support. Seems important not to let actual support
%get up to the boundary of the "known" support, for the positivity to work well.
n=10;
%How many random starts should we attempt:
ntimes=1;
%Set some parameters
ntry=100; %Number of attempts to run the linear prog.
eps=10^(-10); %Terminate when error reaches this size
%Options to control the algorithm
algopts.positivity=0;
algopts.support=1;
%To regularize the support condition
algopts.eps_reg=0.000000;
%This tries to correct for the non-uniqueness that arises due to high order
%of tangency between the torus and the support subspace.
optsexam.k=4;
optsexam.num_disks=5;
%
example = 3;
nplot = 10;
%
for j=1:ntimes
    if example == 1
%A random vector of length N-2n padded with 2n zeros
rng(rand_seed);
x=[zeros(1,n) rand(1,N-2*n) zeros(1,n)];
xx=[x zeros(1,N)];
    elseif example == 2
%A non-random example
x=[abs(sin(4*pi*(1:N-n)/N)),zeros(1,n)];
%Pad to twice the length
xx=[x zeros(1,N)];
%If we create a random collection of "smooth" disks
    elseif example == 3
xsp=(1:2*N)/N;
xx=zeros(1,2*N);
f=zeros(size(xx));
%Set the random number generator
rng(rand_seed);
    for kk=1:optsexam.num_disks
        cc=(rand(3,1)*2-1)*0.7;
        rr=abs((rand*2-1)*0.2);
        R=sqrt((xsp-cc(1)).^2);
        k=optsexam.k;
        xx=xx+(R/rr<1).*((1-(R/rr).^2).^k);
    end;
    xx=xx/max(xx);
    end
%Get the fft
xh=fftn(xx);
S=[eye(N);zeros(N)];
%The support mask for ALG_NEWT
%There is a surprising difference between N and N+1 here
%ones(N+1,1) works better!
algopts.support_mask=[ones(N+1,1); zeros(N-1,1)];
%Compute the intersection dimension of the tang-sp with the supp cond.
    [Tv Nv] = torus_tansp(xx);
    %Here we study the "angle" between the support constraint subspace and the
    %tangent space to the torus.
    [u s v]=svd(S'*Tv);
    d=diag(s); 
    dindex=(d>mincos); %We find the dimension of the space where the 
                       %cosine of the angle is greater than mincos.
    intd=sum(dindex);

magin=abs(xh)';
% figure, plot(magin)
% title(sprintf('Magnitude Fourier data, rand-seed = %d',rand_seed))
%Now  perturb the phases of xh to start our search
%
sigma=sigma_0/sqrt(N);
rng(rand_seed+1);
th=sigma*(rand(1,N-1)-0.5)*2*pi;
figure, hist(2*th/pi)
title('Initial phase errors/\pi/2')
dth=[0 th 0 -fliplr(th)];
%Now we perturb the phases
xh=xh*diag(exp(i*dth));
xh0(j)=xh(1);
xxp=ifftn(xh);
%Plot the initial data
figure
subplot(1,3,1)
plot(xx)
title('Input')
subplot(1,3,2)
plot(abs(fftn(xx)))
title('Input mag FFT') 
xlabel(sprintf('rand-seed = %d',rand_seed))
subplot(1,3,3)
plot(xxp)
title('Init. guess')
%
err=norm(xx-xxp)
x0=xxp';
%ftang=figure; set(ftang,'position',[100,1500,300,300]);
%
tstr=[];
mstr=floor(ntry/10);
D.x=xxp';
for k=1:ntry
    %Newton-like solver
    
    [D,recon,info]=ALG_NEWT(D,magin,algopts);
    x1=D.x;
    err(k)=norm(x1-xx')/norm(xx);
    %if err<eps break
    %end
    %
end
xdif=max(abs(x1-xx'))/max(abs(xx))  
mag_err=norm(abs(fftn(x1))-magin)/norm(magin)
maxx1=max(x1(N+1:2*N))
figure
subplot(1,4,1)
plot(xx)
xlim([1,length(xx)])
title('input')
subplot(1,4,2)
plot(x1)
xlim([1,length(x1)])
title('output')
xlabel(sprintf('R-seed = %d',rand_seed))
subplot(1,4,3)
semilogy(err)
xlabel(sprintf('k = %d',optsexam.k))
title('Recon error')
subplot(1,4,4)
semilogy(abs(abs(fftn(x1'))-abs(fftn(xx))))
xlim([1,length(x1)])
title('Mag error')
xlabel(sprintf('sigma = %g',sigma_0))

end  
