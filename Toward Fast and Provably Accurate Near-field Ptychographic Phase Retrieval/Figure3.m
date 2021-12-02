clear all

%% Code for solving by Wirtinger Flow
ca = cell(1,4);
Tests = 1;

%% Create variables
signalnoiseratio = zeros(1,4);
Error1 = zeros(4,4);
Error2 = zeros(1,15);
d = 30; %size of object
objectX = randn(d,Tests) + 1i*randn(d,Tests);
Maskorg = randn(d,Tests) + 1i*randn(d,Tests);
%Noise = randn(d,d*Tests) + 1i*randn(d,d*Tests);
counter = 0;
for Shifts = 12:2:18
counter = counter + 1;    
ca{counter} = sprintf('K = %d shifts', Shifts);    
K = 0:Shifts-1; %set of shifts
L = 0:d-1;% %set of frequencies
gamma = d/3+1; %Support size of psf in Fourier domain
pointspread = ifft(circshift([ones(gamma,1); zeros(d - gamma,1)],-(gamma-1)/2)); %Creates the low-pass filter pointspread function
Knum = length(K); %Number of shifts
Lnum = length(L); % Number of frequencies
ar = zeros(d,Knum*Lnum);
Ar = zeros((Knum*Lnum)*d,d);
Ar2 = zeros((Knum*Lnum)*d,d);
Y2 = zeros((Knum*Lnum)*d,d);
D = zeros(d,Knum*Lnum);
T = 2000;
s = 0;
for SNR = 20:20:80   
s = s + 1;
Errortest = zeros(1,Tests);
for test = 1:Tests
%Selects the random object and mask for this test
object = objectX(:,test);
maskorg = Maskorg(:,test);



%% Creating measurements

matrixmask = zeros(Lnum,d);
for l = L
   shiftflipp = circshift(reversal(pointspread),-l);
   matrixmask(l+1,:) = conj(shiftflipp(:).*maskorg(:));
end
Yb = zeros(Knum,Lnum);
for k = K
    for l = L
        Yb(k+1,l+1) = abs(innerprod(circshift(conj(matrixmask(l+1,:))',k),object))^2;
    end
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
%noise = Noise(1:Knum,(test-1)*d+1:(test-1)*d+Lnum);
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
[~,~,V] = svd(Z,'econ'); %Computing the largest singular vector
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
z = z - (mu/norm(z0,2)^2)*sum(reshape(temp1(:).*(Matrix*z),d,Knum*Lnum),2);
end
phaseOffset = angle((z'*object)/(object'*object));
objectest = z*exp(1i*phaseOffset);
errorx = 10*log10(norm(objectest - object)^2/norm(object)^2);
[Shifts SNR test T/100 errorx]
Errortest(test) = errorx;

end

signalnoiseratio(s) = SNR;
Error1(s,counter) = mean(rmmissing(Errortest));
end
end

%% Code for solving by Wirtinger Flow

%% Create variables

SNR = 80;
X = 2:2:30;
counter = 0;
for shifts = X
counter = counter + 1;
K = 0:(shifts-1);
Knum = length(K);
ar = zeros(d,Knum*Lnum);
Ar = zeros((Knum*Lnum)*d,d);
Ar2 = zeros((Knum*Lnum)*d,d);
Y2 = zeros((Knum*Lnum)*d,d);
D = zeros(d,Knum*Lnum);



Errortest = zeros(1,Tests);
for test = 1:Tests
object = objectX(:,test);
maskorg = Maskorg(:,test);
%% Creating measurements

matrixmask = zeros(Lnum,d);
for l = 1:Lnum
   shiftflipp = circshift(reversal(pointspread),-l+1);
   matrixmask(l,:) = shiftflipp(:).*maskorg(:);
end
Y = zeros(Knum,Lnum);
for k = K
    for l = L
        Y(k+1,l+1) = abs(innerprod(circshift(conj(matrixmask(l+1,:))',k),object))^2;
    end
end




%noise = Noise(1:Knum,(test-1)*d+1:(test-1)*d+Lnum);
noise = randn(Knum,Lnum) + 1i*randn(Knum,Lnum);
noise = (norm(Y)/10^(SNR/10))*noise/norm(noise);

Y = transpose(Y + noise);
vecY = Y(:);

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
Z = D*((1/(Knum*Lnum))*Matrix.*Y2);
[~,~,V] = svd(Z,'econ');
u = V(:,1);
u = u./abs(u);
z0 = sqrt(lambda)*u;

%% Gradient Descent

mu0 = 0.4;
t0 = 330;
z = z0;
for t = 1:T
temp1 = repmat(transpose(abs(ar'*z).^2 - vecY),d,1);
z = z - (d*min(1 - exp(-t/t0),mu0)/norm(z0,2)^2)*sum(reshape(temp1(:).*(Matrix*z),d,Knum*Lnum),2);
end
phaseOffset = angle((z'*object)/(object'*object));
objectest = z*exp(1i*phaseOffset);
errorx = 10*log10(norm(objectest - object)^2/norm(object)^2);
[shifts SNR test T/100 errorx]
Errortest(test) = errorx;

end
Error2(counter) = mean(rmmissing(Errortest));
end

plot(signalnoiseratio,Error1(:,1),'-b','LineWidth',2)
hold on
plot(signalnoiseratio,Error1(:,2),'--r','LineWidth',2)
hold on
plot(signalnoiseratio,Error1(:,3),':g','LineWidth',2)
hold on
plot(signalnoiseratio,Error1(:,4),'-.k','LineWidth',2)

xlabel({'SNR (in dB)'})
ylabel({'Reconstruction Error (in dB)'})
title({'SNR vs Reconstruction Error'})
xticks(20:10:80) 
legend(ca, 'Location', 'northeast')

figure()

plot(X,Error2,'Marker','o','Color','b','LineWidth',2)
xlabel({'Number of Shifts'})
ylabel({'Reconstruction Error (in dB)'})
title({'Number of Shifts vs Reconstruction Error'})
xticks(2:2:30)  

%% Pre-Assigned Functions


function[f] = innerprod(x,y) %Calculates the inner product of two column vectors of the same size
conjy = conj(y);
f = sum(x(:).*conjy(:));
end
function[f] = reversal(x) %Calculates the reversal of a vector about its first entry
f = circshift(flip(x),1);
end
