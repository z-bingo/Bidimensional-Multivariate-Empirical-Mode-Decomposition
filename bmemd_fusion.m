%This is the demo of fusing multiple images by using bmemd
if exist('im', 'var')
    clear im
end

im1=imread('./IMG/source02_0.tif');
im2=imread('./IMG/source02_1.tif');
color = size(im1,3);

final_img1 = uint8(zeros(size(im1)));
final_img2 = uint8(zeros(size(im1)));

for color_chanel = 1:color
    im(:,:,1)=im1(:,:,color_chanel);
    im(:,:,2)=im2(:,:,color_chanel);
    imf = bmemd(im,8);
    n_imf = size(imf,2);
    [M, N, dim] = size(imf{1,1});
    cor = zeros(M,N,dim);
    var = zeros(M,N,dim);
    fusion_img = zeros(M,N);
    for imf_i=1:n_imf
         if imf_i ~= n_imf
            var = local_var_img(imf{1,imf_i},5);
            for j = 1:dim
                cor(:,:,j) = var(:,:,j) ./ sum(var,3);
            end
        else
            for j = 1:dim
                cor(:,:,j) = imf{1,imf_i}(:,:,j) ./ sum(imf{1,imf_i},3);
            end
        end
        fusion_img = fusion_img + sum(imf{1,imf_i} .* cor,3);
    end
end

imshow(fusion_img, [])

