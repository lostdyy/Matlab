function [ phi ] = Sat( r,sigma)
%UNTITLED3 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

phi = sign(r)*min(abs(r),sigma);

end

