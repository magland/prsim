function [ind]=i2s(sz,i)
%This function converts a linear index into a subscript into the n=length(sz)
%dimensional array with sizes [sz(1),...,sz(n)].
n=length(sz);
%Check to see if index is too large
if i > prod(sz)
    error('index exceeds matrix size')
end
ind=zeros(1,n);
p=prod(sz(1:(n-1)));
for k=1:(n-1)
    in=floor((i-1)/p);
    i=i-in*p;
    ind(n-k+1)=in+1;
    p=p/sz(n-k);
end
ind(1)=i;
end