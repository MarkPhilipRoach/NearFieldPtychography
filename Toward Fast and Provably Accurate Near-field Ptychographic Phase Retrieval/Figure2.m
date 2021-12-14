clear all

%% Code for solving by BlockPR

Tests = 100; %Choose number of tests
%% Assigning variables

d = 102; %Choose the lenghth of the sample
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

delta = 26; %Choose delta
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
maskmatrix = BlockCirculant(masklmatrix,d);

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
noise = randn(Knum,Lnum) + 1i*randn(Knum,Lnum);
noise = (norm(Y)/10^(SNR/10))*noise/norm(noise);
Y = transpose(Y + noise);
%% Constructing y
%Here is the main part of the algorithm, computing the inverse of our
%matrix and multiplying with the vectorization of our noisy measurements

y = maskmatrix\Y(:);
y = y(:);
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
LG = Dmatrix - Xtilde; %Compute the weighted Laplacian
[xrec2, ~, ~] = eigs(LG, 1, 'smallestabs');    % compute smallest eigenvector
xrec2 = xrec2./abs(xrec2); %Normalize this eigenvector
xest = sqrt(diag(X)).*xrec2; %Compute our estimate
phaseOffset = angle( (xest'*object) / (object'*object) ); %Compute the global phase error
xest = xest*exp(1i*phaseOffset); %Fix the global error
errorx = 10*log10(norm(xest - object)^2/norm(object)^2); %Compute the reconstruction error
[delta SNR test errorx toc] %Output results
Errortest(test) = errorx; %Log the reconstruction error
runtime(test) = toc; %End the timer
end
signalnoiseratio(s) = SNR; %Record the SNR used for the test
Error1(s) = mean(Errortest); %Compute the mean of the reconstruction error
runtime1(s) = mean(runtime); %Compute the mean of the runtime
end

%% Code for solving by Wirtinger Flow

%% Dummy variables
signalnoiseratio = zeros(4,1);
Error2 = zeros(4,4);
runtime2 = zeros(4,4);
ar = zeros(d,Knum*Lnum);
Ar = zeros((Knum*Lnum)*d,d);
D = zeros(d,Knum*Lnum);
Errortest = zeros(4,Tests);
runtime = zeros(4,Tests);
%% Iterations
T = [ 250 500 1000 2000]; %Set of Number of Iterations
%T = [100 150 200 250];

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
%% Wirtinger Flow
%Here we now perform our Wirtinger flow gradient descent
mu0 = 0.4; 
t0 = 330;
z = z0;
for t = 1:T(4) %Compute the iterations
arz = ar'*z;
Arz = ar*(arz.*reshape(transpose(abs(arz).^2 - vecY),[],1));
z = z - (min(1 - exp(-t/t0),mu0)/(abs(lambda)))*Arz; %Generate the new iterate
if t==T(1)
phaseOffset = angle((z'*object)/(object'*object)); %Compute the global phase error
objectest = z*exp(1i*phaseOffset); %Fix the global error
errorx = 10*log10(norm(objectest - object)^2/norm(object)^2); %Compute the reconstruction error
[t/100 SNR test  errorx toc] %Output results
Errortest(1,test) = errorx; %Log the reconstruction error
runtime(1,test) = toc; %End the timer
else
end
if t==T(2)
phaseOffset = angle((z'*object)/(object'*object)); %Compute the global phase error
objectest = z*exp(1i*phaseOffset); %Fix the global error
errorx = 10*log10(norm(objectest - object)^2/norm(object)^2); %Compute the reconstruction error
[t/100 SNR test  errorx toc] %Output results
Errortest(2,test) = errorx; %Log the reconstruction error
runtime(2,test) = toc; %End the timer
else
end
if t==T(3)
phaseOffset = angle((z'*object)/(object'*object)); %Compute the global phase error
objectest = z*exp(1i*phaseOffset); %Fix the global error
errorx = 10*log10(norm(objectest - object)^2/norm(object)^2); %Compute the reconstruction error
[t/100 SNR test  errorx toc] %Output results
Errortest(3,test) = errorx; %Log the reconstruction error
runtime(3,test) = toc; %End the timer
else
end
end
phaseOffset = angle((z'*object)/(object'*object)); %Compute the global phase error
objectest = z*exp(1i*phaseOffset); %Fix the global error
errorx = 10*log10(norm(objectest - object)^2/norm(object)^2); %Compute the reconstruction error
[T(4)/100 SNR test errorx toc] %Output results
Errortest(4,test) = errorx; %Log the reconstruction error
runtime(4,test) = toc;%End the timer
end

signalnoiseratio(s) = SNR; %Record the SNR used for the test
for n=1:4
Error2(s,n) = mean(Errortest(n,:)); %Compute the mean of the reconstruction error
runtime2(s,n) = mean(runtime(n,:)); %Compute the mean of the runtime
end
end
%% Plotting figures

%First we plot our first figure, comparing the reconstruction error of the
%two algorithms, applying various numbers of iterations, versus varying levels of SNR

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create multiple lines using matrix input to plot
plot1 = plot(signalnoiseratio,[Error2(:,1) Error2(:,2) Error2(:,3) Error2(:,4) Error1(:,1)],'LineWidth',2,'Parent',axes1);
set(plot1(1),'DisplayName','WF = 250 iters','Color',[0 0 1]);
set(plot1(2),'DisplayName','WF = 500 iters','LineStyle','-.',...
    'Color',[1 0 0]);
set(plot1(3),'DisplayName','WF = 1000 iters','LineStyle',':',...
    'Color',[0 1 0]);
set(plot1(4),'DisplayName','WF = 2000 iters','LineStyle','-.',...
    'Color',[0 0 0]);
set(plot1(5),'DisplayName','BlockPR','Color',[1 0 1]);

% Create ylabel
ylabel({'Reconstruction Error (in dB)'});

% Create xlabel
xlabel({'SNR (in dB)'});

% Create title
title({'SNR vs Reconstruction Error'});

box(axes1,'on');
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'XTick',[20 30 40 50 60 70 80]);
% Create legend
legend(axes1,'show');

%Secondly, we plot the average runtime
%two algorithms, applying various numbers of iterations.


% Create figure
figure2 = figure;

% Create axes
axes1 = axes('Parent',figure2);
hold(axes1,'on');

% Create bar
X = categorical({'BlockPR','WF = 250 Iters','WF = 500 Iters','WF = 1000 Iters', 'WF = 2000 Iters'});
X = reordercats(X,{'BlockPR','WF = 250 Iters','WF = 500 Iters','WF = 1000 Iters', 'WF = 2000 Iters'});
Y = [ mean(runtime1(:,1)) mean(runtime2(:,1)) mean(runtime2(:,2)) mean(runtime2(:,3)) mean(runtime2(:,4))];
bar(X,Y);

% Create ylabel
ylabel('Runtime (in seconds)');

% Create title
title('Iterations vs Runtime');

box(axes1,'on');
hold(axes1,'off');


%% Pre-Assigned Functions

function[f] = reversal(x)
f = circshift(flip(x),1);
end