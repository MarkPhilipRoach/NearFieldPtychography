clear all

Tests = 100; %Choose number of tests
%% Assigning variables

d = 945; %Choose the lenghth of the sample
objectX = randn(d,Tests)+ 1i*randn(d,Tests); %Generate the test samples

%% Dummy variables

%Here we compile the dummy variables that will be used
Errortest = zeros(Tests,1);
runtime = zeros(Tests,1);
signalnoiseratio = zeros(1,4);
Error1 = zeros(4,4);
runtime1 = zeros(1,1);

%% Choose delta
counter = 0;
for delta = [5 14 23 32] %Choose delta
counter = counter + 1; 
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
runtime(test) = toc; %End the timer
phaseOffset = angle( (xest'*object) / (object'*object) ); %Compute the global phase error
xest = xest*exp(1i*phaseOffset); %Fix the global error
errorx = 10*log10(norm(xest - object)^2/norm(object)^2); %Compute the reconstruction error
[delta SNR test errorx toc] %Output results
Errortest(test) = errorx; %Log the reconstruction error

end
signalnoiseratio(s) = SNR; %Record the SNR used for the test
Error1(s,counter) = mean(Errortest); %Compute the mean of the reconstruction error
runtime1(s,counter) = mean(runtime); %Compute the mean of the runtime
end
end

%% Plotting figures

%First we plot our first figure, comparing the reconstruction error with
%the SNR

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create multiple lines using matrix input to plot
plot1 = plot(signalnoiseratio,[Error1(:,1) Error1(:,2) Error1(:,3) Error1(:,4)],'LineWidth',1.5,'Parent',axes1);
set(plot1(1),'DisplayName','\delta = 5','Color',[0 0 1]);
set(plot1(2),'DisplayName','\delta = 14','LineStyle','--','Color',[1 0 0]);
set(plot1(3),'DisplayName','\delta = 23','LineStyle',':','Color',[0 1 0]);
set(plot1(4),'DisplayName','\delta = 32','LineStyle','-.','Color',[0 0 0]);

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

%Secondly, we plot our runtime comparisons versus the delta level

% Create figure
figure2 = figure;

% Create axes
axes1 = axes('Parent',figure2);
hold(axes1,'on');

% Create plot
X = categorical({'5','14','23','32'}); 
X = reordercats(X,{'5','14','23','32'});
Y = [mean(runtime1(:,1)) mean(runtime1(:,2)) mean(runtime1(:,3)) mean(runtime1(:,4))];
plot(X,Y,'Marker','o','LineWidth',2,'Color',[0 0 1]);

% Create ylabel
ylabel({'Runtime (in seconds)'});

% Create xlabel
xlabel({'\delta'});

% Create title
title({'\delta vs Runtime'});

box(axes1,'on');
hold(axes1,'off');

%% Pre-Assigned Functions
function[f] = reversal(x)
f = circshift(flip(x),1);
end
