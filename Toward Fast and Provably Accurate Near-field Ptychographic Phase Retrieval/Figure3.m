clear all

Tests = 100; %Choose number of tests
%% Assigning variables

d = 102; %Choose the lenghth of the sample
objectX = randn(d,Tests)+ 1i*randn(d,Tests); %Generate the test samples

%% Dummy variables

%Here we compile the dummy variables that will be used
Errortest = zeros(4*Tests,2);
runtime = zeros(4*Tests,2);
signalnoiseratio = zeros(1,2);
Error1 = zeros(4,2);
runtime1 = zeros(1,1);
Error2 = zeros(4,2);
runtime2 = zeros(4,2);

ErrortestWF = 10*ones(4*Tests,2);
runtimeWF = zeros(4*Tests,2);
%% Choose delta
delta = 26; %Choose delta
D = d*(2*delta-1);
K = 0:d-1; %set of shifts
L = 0:2*delta-2;% %set of frequencies
Knum = length(K); %Number of shifts
Lnum = length(L); %Number of frequencies
ar = zeros(d,Knum*Lnum);

%% Generating set of point spread functions and masks
%Compute the point spread functions
gamma = floor(d/3+1);
pointspreadP(:,1:Tests) = ifft(circshift([ones(gamma,Tests) ; zeros(d-gamma,Tests)] + 0.01*randn(d,Tests),-floor((gamma-1)/2)));
for t = 1:Tests
    for n=1:d
        pointspreadP(n,Tests+t) = exp(2*pi*1i*randn(1));
    end
end

% Compute the masks

maskorgsmall = randn(delta/2,Tests)+1i*randn(delta/2,Tests);
maskorgM(:,1:Tests) = [maskorgsmall; flip(maskorgsmall); zeros(d-delta,Tests)];
maskorgM(:,Tests+1:2*Tests) = [randn(delta,Tests)+1i*randn(delta,Tests); zeros(d-delta,Tests)];

for counter = 1:2

%% Noise
%Generate the additive noise (not adjusted for magnitude yet) that will be used for the test
Noise = randn(Tests*Knum,Lnum) + 1i*randn(Tests*Knum,Lnum); 

s = 0;
for SNR = 20:20:80 %Adjusting the signal-to-noise ratio
s = s + 1;


%% Starting test
for test = 1:Tests
object = objectX(:,test); %Choosing sample
pointspread = pointspreadP(:,(counter-1)*Tests+test);
maskorg = maskorgM(:,(counter-1)*Tests+test);
stepsize = [2*sqrt(d) 1/d^2]; %Low Pass and Unit Mag
iterations = [5000 1000]; %Low Pass and Unit Mag
%% Constructing measurements
%Construct the convolutional measurements and rearrange the measurement to
%our requried form using Remark 1
Y = zeros(Knum,Lnum);

for k = 1:Knum
    for l = 1:Lnum
    Yrow = cconv(pointspread, circshift(maskorg,k-1).* object,d);
    Y(k,l) = abs(Yrow(mod(k-l,d)+1))^2;
    end
end

%% Adding noise
%Here we add the noise to our measurement for each test
noise = Noise((test-1)*Knum+1:test*Knum,:);
noise = (norm(Y)/10^(SNR/10))*noise/norm(noise);
Y = transpose(Y + noise);

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
runtime(test,counter) = toc; %End the timer
phaseOffset = angle( (xest'*object) / (object'*object) ); %Compute the global phase error
xest = xest*exp(1i*phaseOffset); %Fix the global error
errorx = 10*log10(norm(xest - object)^2/norm(object)^2); %Compute the reconstruction error
[counter SNR test errorx toc] %Output results
Errortest((s-1)*Tests+test,counter) = errorx; %Log the reconstruction error

%% Wirtinger Flow
for n = 1:Knum*Lnum
    ar(:,n) = circshift(conj(matrixmask(L(mod(n-1,Lnum)+1)+1,:))',K(floor((n-1)/Lnum)+1));
end
vecY = reshape(Y,[],1);
mu0 = 0.4; 
t0 = 330;
z0 = xest;
z = z0;
for t=1:iterations(counter)
arz = ar'*z;
Arz = ar*(arz.*reshape(transpose(abs(arz).^2 - vecY),[],1));
z = z - stepsize(counter)*(min(1 - exp(-t/t0),mu0)*Arz/(norm(z0,2)^2)); %Generate the new iterate
phaseOffset = angle((z'*object)/(object'*object)); %Compute the global phase error
objectest = z*exp(1i*phaseOffset); %Fix the global error
errorx = 10*log10(norm(objectest - object)^2/norm(object)^2); %Compute the reconstruction error
ErrortestWF((s-1)*Tests+test,counter) = errorx; %Log the reconstruction error
[counter t/100 SNR test Errortest((s-1)*Tests+test,counter) ErrortestWF((s-1)*Tests+test,counter) toc] %Output results
end
runtimeWF((s-1)*Tests+test,counter) = toc;%End the timer


end
signalnoiseratio(s) = SNR; %Record the SNR used for the test
Error1(s,counter) = mean(rmmissing(Errortest((s-1)*Tests+1:s*Tests,counter))); %Compute the mean of the reconstruction error
runtime1(s,counter) = mean(runtime((s-1)*Tests+1:s*Tests,counter)); %Compute the mean of the runtime
Error2(s,counter) = mean(rmmissing(ErrortestWF((s-1)*Tests+1:s*Tests,counter))); %Compute the mean of the reconstruction error
runtime2(s,counter) = mean(runtimeWF((s-1)*Tests+1:s*Tests,counter)); %Compute the mean of the runtime
end
end


%% Plotting figures
%First we plot our first figure, comparing the reconstruction error with
%the SNR for the low pass filter and symmetric mask

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create multiple lines using matrix input to plot
plot1 = plot(signalnoiseratio,[Error1(:,1) Error2(:,1)],'LineWidth',1.5,'Parent',axes1);
set(plot1(1),'DisplayName','Low Pass PSF, Random Mask ','Color',[0 0 1]);
set(plot1(2),'DisplayName','WF Iteration','LineStyle','--','Color',[1 0 0]);

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

%Secondly, we plot our second figure, comparing the reconstruction error with
%the SNR for the unit magnitude psf and random mask

% Create figure

figure2 = figure;

% Create axes
axes1 = axes('Parent',figure2);
hold(axes1,'on');

% Create multiple lines using matrix input to plot
plot1 = plot(signalnoiseratio,[Error1(:,2) Error2(:,2)],'LineWidth',1.5,'Parent',axes1);
set(plot1(1),'DisplayName','Unit Mag PSF, Symmetric Mask','Color',[0 0 1]);
set(plot1(2),'DisplayName','WF Iteration','LineStyle','--','Color',[1 0 0]);

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

%% Pre-Assigned Functions
function[f] = reversal(x)
f = circshift(flip(x),1);
end