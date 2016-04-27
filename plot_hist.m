figure
subplot(1,2,1)
histogram(intdim1(:,2),20,'BinLimits',[0,max(intdim1(:,2))])
mn1=mean(intdim1(:,2));
title(sprintf('1d, mean = %g',mn1))
subplot(1,2,2)
histogram(intdim2(:,2),20,'BinLimits',[0,max(intdim2(:,2))])
mn2=mean(intdim2(:,2));
title(sprintf('2d, mean = %g',mn2))