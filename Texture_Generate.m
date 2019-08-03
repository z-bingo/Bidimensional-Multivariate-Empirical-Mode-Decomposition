%The code to generate synthetic texture images in paper
clc, close all
f1 = 0.08
f2 = 0.02
n = 0:199;

mask = mat2gray(normpdf(1:200, 100, 80));

s1 = sin(2*pi*f1*n+rand()*pi) + sin(2*pi*f2*n+rand()*pi) + sin(2*pi*f3*n+rand()*pi);
text1 = repmat(s1, 200, 1) + repmat(s1', 1, 200);
text1 = mat2gray(text1);

s2 = sin(2*pi*f2*n+rand()*pi) + sin(2*pi*f1*n+rand()*pi);
text2 = repmat(s2, 200, 1) + repmat(s2', 1, 200);
text2 = mat2gray(text2);
%figure, imshow(text2)
s3 = sin(2*pi*f1*n+rand()*pi) + sin(2*pi*f2*n+rand()*pi);
text3 = repmat(s3, 200, 1) + repmat(s3', 1, 200);
text3 = mat2gray(text3);

mask = repmat(normpdf(1:200, 100, 80), 200, 1)+repmat(normpdf(1:200, 100, 80)', 1, 200);
syn1 = mat2gray(text1);
syn2 = mat2gray(text2+sqrt(0.006)*randn(200,200));
syn3 = mat2gray(text3+sqrt(0.006)*randn(200,200));
% final multi-channel
final = zeros(200, 200, 3);
final(:,:,1) = syn1;
final(:,:,2) = syn2;
final(:,:,3) = syn3;
figure
subplot(2,3,1), imshow(syn1)
subplot(2,3,2), imshow(syn2)
subplot(2,3,3), imshow(syn3)

save('texture_syn.mat', 'final')

