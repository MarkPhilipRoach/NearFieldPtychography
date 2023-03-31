clear all
Tests = 100;


counter = 0;
for SNR = 20:5:40
counter = counter + 1;

colorimage = imread('Ash128.png');
[rows, columns, numberOfColorBands] = size(colorimage);
for s=1:3
tic    
grayimage = colorimage(:, :, s); % Take green channel.
noise = 0*ceil(randn(rows,columns));
for i = 1:rows
    for j = 1:columns
        grayimage(i,j) = grayimage(i,j) + noise(i,j);
    end
end
% Display the original gray scale image.
subplot(2, 3, 1);
imshow(colorimage, []);
fontSize = 15;
title('Original Image', 'FontSize', fontSize, 'Interpreter', 'None');

objectorg = reshape(im2double(grayimage),[],1);

dorg = size(objectorg,1); %Choose the lenghth of the sample


%% Choose delta
delta = 16; %Choose delta
dextend = ceil(dorg/(2*delta-1))*(2*delta-1)-dorg;
%object = [objectorg; objectorg(1:dextend)];
object = [objectorg; ones(dextend,1)];
d = dorg+dextend;
D = d*(2*delta-1);
K = 0:d-1; %set of shifts
L = 0:2*delta-2;% %set of frequencies
Knum = length(K); %Number of shifts
Lnum = length(L); %Number of frequencies
%ar = zeros(d,Knum*Lnum);



%% Generating set of point spread functions and masks
% %Compute the point spread functions
% gamma = floor(d/3+1);
% pointspread = ifft(circshift([ones(gamma,1) ; zeros(d-gamma,1)] + 0.01*randn(d,1),-floor((gamma-1)/2)));
% % Compute the masks
% maskorgsmall = randn(delta/2,1)+1i*randn(delta/2,1);
% maskorg = [maskorgsmall; flip(maskorgsmall); zeros(d-delta,1)];

pointspread = zeros(d,1);
for n = 1:d
    pointspread(n) = exp(-2*pi*1i*n^2/(2*delta-1)); 
end
 %Compute the mask
maskorg = zeros(d,1);
a = max(4,(delta-1)/2);
for n = 1:delta
    maskorg(n) = (exp((-n+1)/a))/((2*delta-1)^(1/4))*exp(2*pi*1i*n^2/(2*delta-1));
end

%% Constructing measurements
%Construct the convolutional measurements and rearrange the measurement to
%our requried form using Remark 1
Y = zeros(Knum,Lnum);

for k = 1:Knum
    %Yrow = cconv(pointspread, circshift(maskorg,k-1).* object,d);
    Yrow = ifft(fft(pointspread).* fft(circshift(maskorg,k-1).* object));
    for l = 1:Lnum     
    Y(k,l) = abs(Yrow(mod(k-l,d)+1))^2;
    end
end


%% Adding noise
%Here we add the noise to our measurement for each test
noise = randn(Knum,Lnum) + 1i*randn(Knum,Lnum); 
noise = (norm(Y)/10^(SNR/10))*noise/norm(noise);
Y = transpose(Y + noise);

%% Constructing phase retrieval masks and matrix Mcheck
% Here we construct a matrix whose rows are the new phase retrieval masks
% generated from our pointspread fucntion and physical mask
matrixmask = zeros(Lnum,d);
for l = 1:Lnum
   shiftflipp = circshift(reversal(pointspread),-l+1);
   matrixmask(l,:) = conj(shiftflipp(:).*maskorg(:));
end

%Here we generate the matrices which will be the blocks of our matrix
%Mcheck
masklmatrix = zeros(2*delta-1,D);
for l = 1:2*delta-1
    for i=1:2*delta-1
        for j=1:delta-l+1
            masklmatrix(i,j+(l-1)*(2*delta-1)) = matrixmask(i,l)*conj(matrixmask(i,j+l-1));
        end
        for j=2*delta-l:2*delta-1
            masklmatrix(i,j+(l-1)*(2*delta-1)) = matrixmask(i,l+1)*conj(matrixmask(i,j+l-2*delta+1));
        end
    end
end

%Construct our block circulant matrix
maskmatrix = BlockCirculant(masklmatrix,d);
%% Constructing y
%Here is the main part of the algorithm, computing the inverse of our
%matrix and multiplying with the vectorization of our noisy measurements

y = maskmatrix\Y(:);
y = y(:);
%% Weighted angular synchronization
%Perform weighted angular synchronization

%Here we form our matrix Xhat
circy = circshift(y,delta-1);
Xhat = zeros(d,d);
for n=1:d
    Xhat(:,n) = circshift([circy((n*(2*delta-1)-(2*delta-2)):n*(2*delta-1)); zeros(d-(2*delta-1),1)],n-delta);
end

%Here we generate our weighted  and degree matrices which will be used to compute our
%weighted Laplacaian
Xhattilde = Xhat.*abs(Xhat);
Weight = abs(Xhat).^2; %Compute weight matrix
Dmatrix = diag(Weight*ones(d,1)); %Compute degree matrix
LG = Dmatrix - Xhattilde; %Compute the weighted Laplacian
[xrec, ~, ~] = eigs(LG, 1, 'smallestabs');    % compute smallest eigenvector
xrec = xrec./abs(xrec); %Normalize this eigenvector
xest = sqrt(diag(Xhat)).*xrec; %Compute our estimate
phaseOffset = angle( (xest'*object) / (object'*object) ); %Compute the global phase error
xest = xest*exp(1i*phaseOffset); %Fix the global error

rgb(:,s) = real(xest(1:dorg));
[SNR s toc]
end

Imageestimate = cat(3, reshape(rgb(:,1),rows,columns), reshape(rgb(:,2),rows,columns), reshape(rgb(:,3),rows,columns));
subplot(2, 3, counter+1);
imshow(Imageestimate, []);
fontSize = 15;
title(sprintf('SNR = %d db',SNR), 'FontSize', fontSize, 'Interpreter', 'None');
%[SNR toc]
end

sgtitle('Image Recovery (128 x 128 pixels, \delta = 16)')
%% Pre-Assigned Functions
function[f] = reversal(x)
f = circshift(flip(x),1);
end

