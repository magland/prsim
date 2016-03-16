function [f,support_mask]=PH_1(opts)

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
f=f+0;
f=f/max(f(:));

end