clear all

Tests = 1; %Choose number of tests
ca = cell(1,5); %Generate blank legend for figure
%% Assigning variables

d = 30; %Choose the lenghth of the sample
objectX = rand(d,Tests) + 1i*randn(d,Tests); %Generate the test samples


%% Dummy variables

%Here we compile the dummy variables that will be used
Errortest = zeros(Tests,1);
runtime = zeros(Tests,1);
Xtilde = zeros(d,d);
Weight = zeros(d,d);
Dmatrix = zeros(d,d);
signalnoiseratio = zeros(4,1);
Error1 = zeros(4,1);
runtime1 = zeros(1,1);
%% Choose delta

delta = 8;%Choose delta
ca{5} = sprintf('%d BlockPR', delta); 
D = d*(2*delta-1);
K = 0:d-1; %set of shifts
L = 0:2*delta-2;% %set of frequencies
Knum = length(K); %Number of shifts
Lnum = length(L); %Number of frequencies
%% Generating point spread function and mask
%Compute the point spread function
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

%% Noise
%Generate the additive noise (not adjusted for magnitude yet) that will be used for the test
Noise = randn(Tests*Knum,Lnum) + 1i*randn(Tests*Knum,Lnum); 

s = 0;
for SNR = 20:20:80 %Adjusting the signal-to-noise ratio
s = s + 1;
%% Starting test
for test = 1:Tests
object = objectX(:,test); %Choosing sample
tic %Start the timer

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

%Construct our matrix
maskmatrix = repmat(masklmatrix,Knum,1);
for i=1:Knum
    maskmatrix((i-1)*(2*delta-1)+1:i*(2*delta-1),:) = circshift(maskmatrix((i-1)*(2*delta-1)+1:i*(2*delta-1),:) ,(i-1)*(2*delta-1),2);
end

%% Constructing measurements
%Construct the convolutional measurements and rearrange the measurement to
%our requried form
Yconv = zeros(Knum,Lnum);
Y = zeros(Knum,Lnum);
for k = 1:d
    Yrow = cconv(pointspread, circshift(maskorg,-k+1).* object,d);
    for l = 1:d
        Yconv(k,l) = abs(Yrow(l))^2;
    end
end
for k = 1:Knum
    for l = 1:Lnum
        Y(k,l) = Yconv(mod(-k+1,d)+1,mod(-l+k,2*delta-1)+1);
    end
end

%% Adding noise
%Here we add the noise to our measurement for each test
noise = Noise((test-1)*Knum+1:test*Knum,:);
noise = (norm(Y)/10^(SNR/10))*noise/norm(noise);
Y = transpose(Y + noise);
%% Constructing y
%Here is the main part of the algorithm, computing the inverse of our
%matrix and multiplying with the vectorization of our noisy measurements

y = maskmatrix\Y(:);

%% Weighted angular synchronization
%Perform weighted angular synchronization

%Here we form our matrix Xhat
circy = circshift(y,delta-1);
X = zeros(d,d);
for n=1:d
    X(:,n) = circshift([circy((n*(2*delta-1)-(2*delta-2)):n*(2*delta-1)); zeros(d-(2*delta-1),1)],n-delta);
end

%Here we generate our weighted  and degree matrices which will be used to compute our
%weighted Laplacaian
for i=1:d
    for j=1:d
        if X(i,j)==0
        else
            if i==j
            else
            Xtilde(i,j) = X(i,j)/norm(X(i,j))*abs(X(i,j))^2;
            Weight(i,j) = abs(X(i,j))^2;
            end
        end
    end
    Dmatrix(i,i) = sum(Weight(i,:));
end
L = Dmatrix - Xtilde; %Compute the weighted Laplacian
[xrec2, ~, ~] = eigs(L, 1, 'smallestabs');    % compute smallest eigenvector
xrec2 = xrec2./abs(xrec2); %Normalize this eigenvector
xest = sqrt(diag(X)).*xrec2; %Compute our estimate
phaseOffset = angle( (xest'*object) / (object'*object) ); %Compute the global phase error
xest = xest*exp(1i*phaseOffset); %Fix the global error
errorx = 10*log10(norm(xest - object)^2/norm(object)^2); %Compute the reconstruction error
[delta SNR test errorx toc] %Output results
Errortest(test) = errorx; %Log the reconstruction error
runtime(test) = toc; %End the timer
end
signalnoiseratio(s) = SNR;
Error1(s) = mean(Errortest);
runtime1(s) = mean(runtime);
end
runtime = zeros(1,Tests);
%% Code for solving by Wirtinger Flow

%% Create variables

signalnoiseratio = zeros(4,1);
Error2 = zeros(4,5);
runtime2 = zeros(4,5);
nc = 0;
 
pointspread = zeros(d,1);
for n = 1:d
    pointspread(n) = exp(2*pi*1i*n^2/(2*delta-1));
end
maskorg = zeros(d,1);
for n = 1:delta
    maskorg(n) = (exp((-n+1)/a))/((2*delta-1)^(1/4))*exp(-2*pi*1i*n^2/(2*delta-1));
end
K = 0:d-1; %set of shifts
L = 0:2*delta-2;% %set of frequencies
Knum = length(K);
Lnum = length(L);
Noise = randn(Tests*Knum,Lnum) + 1i*randn(Tests*Knum,Lnum); 
ar = zeros(d,Knum*Lnum);
Ar = zeros((Knum*Lnum)*d,d);
Ar2 = zeros((Knum*Lnum)*d,d);
Y2 = zeros((Knum*Lnum)*d,d);
D = zeros(d,Knum*Lnum);
for T = [250 500 1000 2000]
nc = nc + 1; 
ca{nc} = sprintf('%d iters', T);    
s = 0;
for SNR = 20:20:80   
s = s + 1;
Errortest = zeros(1,Tests);
runtime = zeros(1,Tests);
for test = 1:Tests
object = objectX(:,test);

%% Creating measurements
matrixmask = zeros(Lnum,d);
for l = L
   shiftflipp = circshift(reversal(pointspread),-l);
   matrixmask(l+1,:) = conj(shiftflipp(:).*maskorg(:));
end
Yconv = zeros(Knum,d);
Y = zeros(Knum,Lnum);
for k = 1:d
    Yrow = cconv(pointspread, circshift(maskorg,-k+1).* object,d);
    for l = 1:d
        Yconv(k,l) = abs(Yrow(l))^2;
    end
end
for k = K
    for l = 1:Lnum
        Y(k+1,l) = Yconv(mod(-k,d)+1,mod(-l+k+1,Lnum)+1);
    end
end
%% Adding Noise

noise = randn(Knum,Lnum) + 1i*randn(Knum,Lnum);
noise = (norm(Y)/10^(SNR/10))*noise/norm(noise);

Y = transpose(Y + noise);
vecY = Y(:);

tic
for n = 1:Knum*Lnum
    ar(:,n) = circshift(conj(matrixmask(L(mod(n-1,Lnum)+1)+1,:))',K(floor((n-1)/Lnum)+1));
        Ar((n-1)*d+1:n*d,:) = repmat(ar(:,n),1,d);
        Ar2((n-1)*d+1:n*d,:) = repmat(ar(:,n)',d,1);
        Y2((n-1)*d+1:n*d,1:d) = vecY(n);
        for i = 1:d
            D(i,(n-1)*d+i) = mod(ceil((i+n-2)/d)+1 ,2);
        end
end
Matrix = (Ar.*Ar2);

%% Initialization

lambda = d*sum(vecY)/sum(vecnorm(ar));
Z = D*(Matrix.*Y2)*(1/(Knum*Lnum));
[~,~,V] = svd(Z,'econ');
u = V(:,1);
u = u./abs(u);
z0 = sqrt(lambda)*u;
%% Gradient Descent

mu0 = 0.4;
t0 = 330;
z = z0;
for t = 1:T
mu = d*min(1 - exp(-t/t0),mu0);
temp1 = repmat(transpose(abs(ar'*z).^2 - vecY),d,1);
z = z - (mu/(d*abs(lambda)))*sum(reshape(temp1(:).*(Matrix*z),d,Knum*Lnum),2);
end
phaseOffset = angle((z'*object)/(object'*object));
objectest = z*exp(1i*phaseOffset);
errorx = 10*log10(norm(objectest - object)^2/norm(object)^2);
[T/100 SNR test (t/T)*100 errorx toc]
Errortest(test) = errorx;
runtime(test) = toc;

end

signalnoiseratio(s) = SNR;
Error2(s,nc) = mean(Errortest);
runtime2(s,nc) = mean(runtime);
end
end

plot(signalnoiseratio,Error2(:,1),'-b','LineWidth',2)
hold on
plot(signalnoiseratio,Error2(:,2),'-.r','LineWidth',2)
hold on
plot(signalnoiseratio,Error2(:,3),':g','LineWidth',2)
hold on
plot(signalnoiseratio,Error2(:,4),'-.k','LineWidth',2)
hold on
plot(signalnoiseratio,Error1(:,1),'m','LineWidth',2)
hold on

xlabel({'SNR (in dB)'})
ylabel({'Reconstruction Error (in dB)'})
title({'SNR vs Reconstruction Error'})
xticks(20:10:80) 
legend(ca, 'Location', 'northeast')

figure()

X = categorical({'BlockPR','WF = 250 Iters','WF = 500 Iters','WF = 1000 Iters', 'WF = 2000 Iters'});
X = reordercats(X,{'BlockPR','WF = 250 Iters','WF = 500 Iters','WF = 1000 Iters', 'WF = 2000 Iters'});
Y = [ mean(runtime1(:,1)) mean(runtime2(:,1)) mean(runtime2(:,2)) mean(runtime2(:,3)) mean(runtime2(:,4))];
bar(X,Y)
title('Iterations vs Runtime')
ylabel('Runtime (in seconds)')


%% Pre-Assigned Functions

function[f] = reversal(x)
f = circshift(flip(x),1);
end
