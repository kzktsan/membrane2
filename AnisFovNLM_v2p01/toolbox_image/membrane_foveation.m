D = dir('/Users/satokazuki/Desktop/zikken4/zikken4_cellmem/*png')
addpath('/Users/satokazuki/Desktop/zikken4/zikken4_cellmem');
n = 256;
%m = 15; %width of filters
%p = 15; % number of filters
m = 41;
p = 10;
sigma = linspace(0.05,10,p);
H = zeros(m,m,p);
for i=1:p
    H(:,:,i) = compute_gaussian_filter([m m],sigma(i)/n,[n n]);
end
x = linspace(-1,1,n);
[Y,X] = meshgrid(x,x);
R = sqrt(X.^2 + Y.^2);
I = round(rescale(R,1,p));
for k=1:length(D)
    file_name = D(k).name;
    image_name = strsplit(file_name, '.');
    M = load_image(char(image_name(1)),n);
    %imageplot(M);
    M1 = perform_adaptive_filtering(M,H,I);
    M2 = uint8(M1);
    save_name = strcat('/Users/satokazuki/Desktop/zikken4fov/zikken4fov_cellmem/',image_name(1),'.png');
    imwrite(M2, char(save_name));
end