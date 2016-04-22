function [ cind ] = conjind(sz,ind )
%Finds the conjugate index to the location i2s(sz,ind) in the array with
%dimensions sz
n=length(sz);
%Unflatten the index
sind=i2s(sz,ind);
%Find the conjugate relative to sz
scind=modm(2*ones(1,n)-sind,sz);
%Flatten the conjugate index
cind=s2i(sz,scind);
end

