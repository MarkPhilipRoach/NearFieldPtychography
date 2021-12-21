clear all

%% Code for solving by Wirtinger Flow

Tests = 100; %Choose number of tests

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
ca = cell(1,4);
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

tic %Start timer
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
end
lambda = d*sum(vecY)/sum(vecnorm(ar));
Z = zeros(d,d);
for n = 1:Knum*Lnum
    Z = Z + vecY(n)*ar(:,n)*ar(:,n)'*1/(Knum*Lnum);
end
[u, ~, ~] = eigs(Z, 1, 'largestabs'); %Compute largest eigenvector
z0 = sqrt(lambda)*u/norm(u); %Set initial estimate

%% Gradient Descent

mu0 = 0.4;
t0 = 330;
z = z0;
for t = 1:T %Compute the iterations
arz = ar'*z;
Arz = ar*(arz.*reshape(transpose(abs(arz).^2 - vecY),[],1));
z = z - (min(1 - exp(-t/t0),mu0)/(abs(lambda)))*Arz; %Generate the new iterate
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
end

lambda = d*sum(vecY)/sum(vecnorm(ar));
Z = zeros(d,d);
for n = 1:Knum*Lnum
    Z = Z + vecY(n)*ar(:,n)*ar(:,n)'*1/(Knum*Lnum);
end
[u, ~, ~] = eigs(Z, 1, 'largestabs'); %Compute largest eigenvector
z0 = sqrt(lambda)*u/norm(u); %Set initial estimate

%% Gradient Descent

mu0 = 0.4;
t0 = 330;
z = z0;
for t = 1:T %Compute the iterations
arz = ar'*z;
Arz = ar*(arz.*reshape(transpose(abs(arz).^2 - vecY),[],1));
z = z - (min(1 - exp(-t/t0),mu0)/(abs(lambda)))*Arz; %Generate the new iterate
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

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create multiple lines using matrix input to plot
plot1 = plot(signalnoiseratio,Error1,'LineWidth',2,'Parent',axes1);
set(plot1(1),'DisplayName','K = 30 shifts','Color',[0 0 1]);
set(plot1(2),'DisplayName','K = 45 shifts','LineStyle','--','Color',[1 0 0]);
set(plot1(3),'DisplayName','K = 60 shifts','LineStyle',':','Color',[0 1 0]);
set(plot1(4),'DisplayName','K = 75 shifts','LineStyle','-.','Color',[0 0 0]);

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

%Secondly, we plot the reconstruction error for fixed SNR, amongst varying
%number of shifts taken

% Create figure
figure2 = figure;

% Create axes
axes1 = axes('Parent',figure2);
hold(axes1,'on');

% Create plot
plot(X,Error2,'Marker','o','LineWidth',2,'Color',[0 0 1]);

% Create ylabel
ylabel({'Reconstruction Error (in dB)'});

% Create xlabel
xlabel({'Number of Shifts'});

% Create title
title({'Number of Shifts vs Reconstruction Error'});

box(axes1,'on');
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'XTick',[15 20 25 30 35 40 45 50 55 60 65 70 75 80 85]);

%% Pre-Assigned Functions

function[f] = reversal(x) %Calculates the reversal of a vector about its first entry
f = circshift(flip(x),1);
end
