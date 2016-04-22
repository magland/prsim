function [ ind ]=s2i(sz,sub)
%This function converts a  linear index into subscripts for an array of
%size sz
%
n=length(sz);
%Check to see if subscripts are allowable
if ~(prod(sub <= sz))
    error('Subscripts are too large')
end
ind=0;
for j=1:(n-1)  
p=prod(sz(1:(n-j)));
ind=ind+p*(sub(n+1-j)-1); 
end
ind=ind+sub(1);
end