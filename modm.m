function y=modm(x,m)
%This function converts a list of integers x into the same list mod m.
%m can either be a single integer, or a list the same length as x.
%
%It differs from the usual function in that mod(m,m)=m (rather than 0)
%
y=mod(x,m);
%Change 0s to ms
y=y+m.*( y == 0);
end