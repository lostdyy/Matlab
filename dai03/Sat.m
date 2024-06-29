function [ phi ] = Sat( r,sigma)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明

phi = sign(r)*min(abs(r),sigma);

end

