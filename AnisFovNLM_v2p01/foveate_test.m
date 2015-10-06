n = 256;
%M = load_image('membrane84', n);
M = load_image('lena', n);
%imageplot(M)
m = 41; % width of filers
p = 10; % number of filters
sigma = linspace(0.05,10,p);
H = zeros(m,m,p);
for i=1:p
    H(:,:,i) = compute_gaussian_filter([m m],sigma(i)/n,[n n]);
end
x = linspace(-1,1,n);
[Y,X] = meshgrid(x,x);
R = sqrt(X.^2 + Y.^2);
I = round(rescale(R,1,p));
M1 = perform_adaptive_filtering(M,H,I);
%imageplot({M M1},{'Original' 'Foveated'});
%imageplot(M1);
M2 = uint8(M1);
imwrite(M2, '/Users/satokazuki/Desktop/testest4.png');