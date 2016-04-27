function [f,support_mask]=PH_1(opts)
%Initialize the random number generator, if needed
if (isfield(opts,'seed')) rng(opts.seed); end;
%
oversamp=opts.oversamp;
N1=opts.N; N2=opts.N;
k=opts.k;

N1b=N1*oversamp;
N2b=N2*oversamp;
[xx,yy]=ndgrid(linspace(-oversamp,oversamp,N1b),linspace(-oversamp,oversamp,N2b));
support_mask=(abs(xx)<=1).*(abs(yy)<=1);

f=zeros(size(xx));
for kk=1:opts.num_disks
    cc=(rand(2,1)*2-1)*0.7;
    rr=abs((rand*2-1)*0.2);
    R=sqrt((xx-cc(1)).^2+(yy-cc(2)).^2);
    f=f+(R/rr<1).*(1-(R/rr).^2).^k;
end;
%Compute more accurate support mask
if ((isfield(opts,'xsupp_msk')) & ~(opts.xsupp_msk == 0))
    suppx=(f> 10^(-12));
    for j=1:opts.xsupp_msk
        suppx=twod_nbhd(suppx);
    end
    support_mask=suppx.*support_mask;
end  
f=f+0;
f=f/max(f(:));

end