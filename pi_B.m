function f2=pi_B(f0,opts)
if opts.positivity
    f0=f0.*(f0>=0);
end;
if opts.support
    f0=f0.*opts.support_mask;
end;
if ~(isfield(opts,'eps_reg')) opts.eps_reg=0; end;
f2=f0/(1+opts.eps_reg);
end