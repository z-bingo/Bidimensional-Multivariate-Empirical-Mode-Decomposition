function [ var ] = local_var_img( in_img, Ni )
%计算灰度图像的局部方差
%输入参数：in_img 要计算局部方差的灰度图像
%         Ni 若为某一固定的常数 表示全局计算局部方差的范围相同
%            也可为与图像等大小的矩阵  表示计算每一点的局部方差范围
%输出参数：var 与图像等大小  表示每一点的局部方差

%先初始化局部方差矩阵  全为0
var = zeros(size(in_img));
%如果Ni是一个数字  那么就讲它的size改为和图像大小相同的常数矩阵
%in_img = im2uint8(in_img);
%计算
mean_img = imfilter(in_img,ones(Ni,Ni)/(Ni*Ni),'replicate');
var = (in_img - mean_img) .^ 2;
end

