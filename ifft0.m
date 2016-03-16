function f=ifft0(fhat)
f=real(fftshift(ifftn(fftshift(fhat))));
end