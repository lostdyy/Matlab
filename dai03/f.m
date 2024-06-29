function [ fx ] = f( x )
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
  fx(1,:) = 0.3*sin(x(1));
  fx(2,:) = 0.2*sin(x(2));

end

