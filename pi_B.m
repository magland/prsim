function f2=pi_B(f0,opts)
if opts.positivity
    f0=f0.*(f0>=0);
end;
if opts.support
    f0=f0.*opts.support_mask;
end;
f2=f0;
end