function resid=compute_residual(f,ref)
resid=sqrt(sum((f(:)-ref(:)).^2))/sqrt(sum(ref(:).^2));
end
