function [ var ] = local_var_img( in_img, Ni )
%����Ҷ�ͼ��ľֲ�����
%���������in_img Ҫ����ֲ�����ĻҶ�ͼ��
%         Ni ��Ϊĳһ�̶��ĳ��� ��ʾȫ�ּ���ֲ�����ķ�Χ��ͬ
%            Ҳ��Ϊ��ͼ��ȴ�С�ľ���  ��ʾ����ÿһ��ľֲ����Χ
%���������var ��ͼ��ȴ�С  ��ʾÿһ��ľֲ�����

%�ȳ�ʼ���ֲ��������  ȫΪ0
var = zeros(size(in_img));
%���Ni��һ������  ��ô�ͽ�����size��Ϊ��ͼ���С��ͬ�ĳ�������
%in_img = im2uint8(in_img);
%����
mean_img = imfilter(in_img,ones(Ni,Ni)/(Ni*Ni),'replicate');
var = (in_img - mean_img) .^ 2;
end

