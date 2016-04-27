function [ supp_plus ] = twod_nbhd( supp0 )
%This takes a 2d support mask and creates the support mask increased by 1
%pixel in each direction
supp1=supp0;
supp1=supp1+circshift(supp0,1,1);
supp1=supp1+circshift(supp1,1,2);
supp1=supp1+circshift(supp1,-1,1);
supp1=supp1+circshift(supp1,-1,2);
supp_plus=(supp1 >0);
end

