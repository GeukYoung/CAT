function g = Imfilt(f,D0,type,n)
% f : input
% D0 : cutoff
% type : 1 - ideal LPF
%        2 - butter LPF
%        3 - gaussian HPF
%        4 - ideal HPF
%        5 - butter HPF
%        6 - gaussian HPF

if nargin < 4
    n = 2;
end

[M,N]=size(f);
F = fft2(double(f));
u = 0:(M-1);
v = 0:(N-1);
idx = find(u>M/2);
u(idx)=u(idx)-M;
idy = find(v>N/2);
v(idy) = v(idy)-N;
[V,U] = meshgrid(v,u);
D = sqrt(U.^2 + V.^2);

switch type
    case 1 % ideal LPF
        H = double(D<=D0);
    case 2 % butter LPF
        H = 1./(1 + (D./D0).^(2*n));
    case 3 % gaussian HPF
        H = exp(-(D.^2)./(2*(D0^2)));
    case 4 % ideal HPF
        H = double(D<=D0);
        H = 1 - H;
    case 5 % butter HPF
        H = 1./(1 + (D0./D).^(2*n));
    case 6 % gaussian HPF
        H = 1 - exp(-(D.^2)./(2*(D0^2)));
end

G = H.*F;
g = real(ifft2(double(G)));
% subplot(1,2,1); imshow(f); title('input');
% subplot(1,2,2); imshow(g,[]); title('output');
end