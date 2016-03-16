function f2=pi_A(f0,opts)
f1=fft0(f0);
f1=opts.u.*exp(i*angle(f1));
f1=real(ifft0(f1));
%f2=f0+opts.alpha_A*(f1-f0);
f2=f1;
end