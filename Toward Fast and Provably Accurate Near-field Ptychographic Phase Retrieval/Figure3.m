clear all

%% Code for solving by Wirtinger Flow

Tests = 1; %Choose number of tests
ca = cell(1,4);
%% Assigning variables
d =  102; %size of object
objectX = randn(d,Tests) + 1i*randn(d,Tests);  %Generate the test samples
Maskorg = randn(d,Tests) + 1i*randn(d,Tests);  %Generate the test masks
T = 2000; %Set number of iterations
%% Dummy variables
runtime = zeros(1,Tests);
signalnoiseratio = zeros(1,4);
Error1 = zeros(4,4);
Error2 = zeros(1,15);
runtime2 = zeros(4,5);
Errortest = zeros(1,Tests);
%% Shifts
counter = 0;
for Shifts = 30:15:75 %Set number of shifts
counter = counter + 1;    
ca{counter} = sprintf('K = %d shifts', Shifts); %Update legend  
K = 0:Shifts-1; %set of shifts
L = 0:d-1;% %set of frequencies
Knum = length(K); %Number of shifts
Lnum = length(L); % Number of frequencies

%% Dummy variables part deux
ar = zeros(d,Knum*Lnum);
Ar = zeros((Knum*Lnum)*d,d);
D = zeros(d,Knum*Lnum);

%% Noise
s = 0;
for SNR = 20:20:80   
s = s + 1;

%% Starting test
for test = 1:Tests
object = objectX(:,test); %Choosing sample
%% Generating point spread function and mask
%Compute the point spread function
gamma = d/3+1; %Support size of psf in Fourier domain
pointspread = ifft(circshift([ones(gamma,1); zeros(d - gamma,1)],-(gamma-1)/2)); %Creates the low-pass filter point-spread function
 %Compute the mask
maskorg = Maskorg(:,test); %Choosing mask


tic 
%% Constructing measurements
%Construct the convolutional measurements and rearrange the measurement to
%our requried form
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

%% Adding noise
%Here we add the noise to our measurement for each test
noise = randn(Knum,Lnum) + 1i*randn(Knum,Lnum);
noise = (norm(Y)/10^(SNR/10))*noise/norm(noise);
vecY = reshape(transpose(Y + noise),[],1);
%% Initialization
for n = 1:Knum*Lnum
    ar(:,n) = circshift(conj(matrixmask(L(mod(n-1,Lnum)+1)+1,:))',K(floor((n-1)/Lnum)+1));
    Ar((n-1)*d+1:n*d,:) = ar(:,n)*ar(:,n)';
    for i = 1:d
            D(i,(n-1)*d+i) = mod(ceil((i+n-2)/d)+1 ,2);
    end
end

lambda = d*sum(vecY)/sum(vecnorm(ar));
Z = D*(Ar.*reshape(repmat(reshape(transpose(repmat(vecY,1,d)),1,[]),1,d),[],d))*(1/(Knum*Lnum));

[u, ~, ~] = eigs(Z, 1, 'largestabs','Tolerance',1e-4); %Compute largest eigenvector
z0 = sqrt(lambda)*u/norm(u); %Set initial estimate

%% Gradient Descent

mu0 = 0.4;
t0 = 330;
z = z0;
for t = 1:T %Compute the iterations
if mod(t-1,50)==0
Arz = Ar*z;
else 
end
z = z - (min(1 - exp(-t/t0),mu0)/(abs(lambda)))*sum(reshape(reshape(repmat(transpose(abs(ar'*z).^2 - vecY),d,1),[],1).*Arz,d,Knum*Lnum),2);  %Generate the new iterate
end
phaseOffset = angle((z'*object)/(object'*object)); %Compute the global phase error
objectest = z*exp(1i*phaseOffset); %Fix the global error
errorx = 10*log10(norm(objectest - object)^2/norm(object)^2); %Compute the reconstruction error
[Shifts SNR test T/100 errorx toc] %Output results
Errortest(test) = errorx; %Log the reconstruction error

end

signalnoiseratio(s) = SNR; %Record the SNR used for the test
Error1(s,counter) = mean(Errortest); %Compute the mean of the reconstruction error
end
end

%% Code for solving by Wirtinger Flow

%% Create variables

SNR = 80; %Fix SNR
counter = 0;
X = 15:5:85;
for Shifts = X %Determine number of shifts
counter = counter + 1;
K = 0:Shifts-1; %Set of shifts
Knum = length(K); %Number of shifts

%% Dummy variables
ar = zeros(d,Knum*Lnum);
Ar = zeros((Knum*Lnum)*d,d);
Ar2 = zeros((Knum*Lnum)*d,d);
Y2 = zeros((Knum*Lnum)*d,d);
D = zeros(d,Knum*Lnum);
Errortest = zeros(1,Tests);


%% Starting test
for test = 1:Tests
object = objectX(:,test); %Choosing sample
maskorg = Maskorg(:,test); %Choosing mask
tic
%% Constructing measurements
%Construct the convolutional measurements and rearrange the measurement to
%our requried form
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

%% Adding noise
%Here we add the noise to our measurement for each test
noise = randn(Knum,Lnum) + 1i*randn(Knum,Lnum);
noise = (norm(Y)/10^(SNR/10))*noise/norm(noise);
Y = transpose(Y + noise);
vecY = Y(:);
%% Initialization

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

lambda = d*sum(vecY)/sum(vecnorm(ar));
Z = D*(1/(Knum*Lnum)*Matrix.*Y2);
[~,~,V] = svd(Z,'econ'); %Compute the SVD of Z
u = V(:,1);%Compute the largest eigenvector
u = u/norm(u); %Normalize eigenvector
z0 = sqrt(lambda)*u; %Generate the initialization

%% Gradient Descent

mu0 = 0.4;
t0 = 330;
z = z0;
for t = 1:T %Compute the iterations
mu = min(1 - exp(-t/t0),mu0); %Compute the stepsize
temp1 = repmat(transpose(abs(ar'*z).^2 - vecY),d,1);
z = z - (mu/abs(lambda))*sum(reshape(temp1(:).*(Matrix*z),d,Knum*Lnum),2); %Generate the new iterate
end
phaseOffset = angle((z'*object)/(object'*object)); %Compute the global phase error
objectest = z*exp(1i*phaseOffset); %Fix the global error
errorx = 10*log10(norm(objectest - object)^2/norm(object)^2); %Compute the reconstruction error
[Shifts SNR test T/100 errorx toc] %Output results
Errortest(test) = errorx; %Log the reconstruction error

end

Error2(counter) = mean(Errortest); %Compute the mean of the reconstruction error
end
%% Plotting figures

%First we plot our first figure, comparing the reconstruction, applying various numbers of shifts, versus varying levels of SNR
plot(signalnoiseratio,Error1(:,1),'-b','LineWidth',2)
hold on
plot(signalnoiseratio,Error1(:,2),'--r','LineWidth',2)
hold on
plot(signalnoiseratio,Error1(:,3),':g','LineWidth',2)
hold on
plot(signalnoiseratio,Error1(:,4),'-.k','LineWidth',2)

xlabel({'SNR (in dB)'}) %Generate label for x-axis
ylabel({'Reconstruction Error (in dB)'}) %Generate label for y-axis
title({'SNR vs Reconstruction Error'}) %Generate title
xticks(20:10:80) 
legend(ca, 'Location', 'northeast') %Generate the legend

figure() %Start new figure

%Secondly, we plot the reconstruction error for fixed SNR, amongst varying
%number of shifts taken

plot(X,Error2,'Marker','o','Color','b','LineWidth',2)
xlabel({'Number of Shifts'})
ylabel({'Reconstruction Error (in dB)'})
title({'Number of Shifts vs Reconstruction Error'})
xticks(X)  

%% Pre-Assigned Functions

function[f] = reversal(x) %Calculates the reversal of a vector about its first entry
f = circshift(flip(x),1);
end
