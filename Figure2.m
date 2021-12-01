clear all

Tests = 1;
ca = cell(1,5);
%% Assigning variables

d = 30;
signalnoiseratio = zeros(4,1);
Error1 = zeros(4,1);
runtime1 = zeros(4,1);
nc = 0;
objectX = rand(d,Tests) + 1i*randn(d,Tests);
delta = 8;
ca{5} = sprintf('%d BlockPR', delta); 
nc = nc + 1; 
D = d*(2*delta-1);
a = max(4,(delta-1)/2);
Noise = randn(Tests*d,2*delta-1) + 1i*randn(Tests*d,2*delta-1);
%% Choice of point spread function & mask


pointspread = zeros(d,1);
for n = 1:d
    pointspread(n) = exp(2*pi*1i*n^2/(2*delta-1));
end
maskorg = zeros(d,1);
for n = 1:delta
    maskorg(n) = (exp((-n+1)/a))/((2*delta-1)^(1/4))*exp(-2*pi*1i*n^2/(2*delta-1));
end

%% Constructing block measurment matrix M'
s = 0;
for SNR = 20:20:80
s = s + 1;
Errortest = zeros(1,Tests);
runtime = zeros(1,Tests);
for test = 1:Tests
x = objectX(:,test);
conjx = conj(x);    
tic

matrixmask = zeros(2*delta-1,d);
for l = 1:2*delta-1
   shiftflipp = circshift(reversal(pointspread),-l);
   matrixmask(l,:) = conj(shiftflipp(:).*maskorg(:));
end

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
maskmatrix = repmat(masklmatrix,d,1);
for i=1:d
    maskmatrix((i-1)*(2*delta-1)+1:i*(2*delta-1),:) = circshift(maskmatrix((i-1)*(2*delta-1)+1:i*(2*delta-1),:) ,(i-1)*(2*delta-1),2);
end



%% Constructing b

b = zeros(d,2*delta-1);
for i = 1:d
    for l = 1:2*delta-1
         b(i,l) = abs(innerprod(matrixmask(l,:),circshift(x,-i+1)))^2;
    end
end

%% Adding noise

noise = Noise((test-1)*d+1:test*d,:);
noise = (norm(b)/10^(SNR/10))*noise/norm(noise);
b = b' + noise';
%% Constructing y

y = maskmatrix\b(:);


%% Angular Synchronization
Xtrue = x*x';
circy = circshift(y,delta-1);
X = zeros(d,d);
for n=1:d
    X(:,n) = [circy((n*(2*delta-1)-(2*delta-2)):n*(2*delta-1)); zeros(d-(2*delta-1),1)];
    X(:,n) = circshift(X(:,n),n-delta);
end
% errorX = norm(X - Xtrue)^2/norm(Xtrue)^2
% mags = sqrt( diag(X) );
% [xrec, ~, ~] = eigs(X, 1, 'LM');    % compute leading eigenvector
% xrec = xrec./abs(xrec);
% xest = sqrt(diag(X)) .* xrec;
% phaseOffset = angle( (xest'*x) / (x'*x) );
% xest = xest*exp(1i*phaseOffset);
% [x xest];
% errorx = 10*log10(norm(xest - x)^2/norm(x)^2);

Xtilde = zeros(d,d);
Weight = zeros(d,d);
Dmatrix = zeros(d,d);
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
L = Dmatrix - Xtilde;
[xrec2, ~, ~] = eigs(L, 1, 'smallestabs');    % compute smallest eigenvector
xrec2 = xrec2./abs(xrec2);
xest = sqrt(diag(X)).*xrec2;
phaseOffset = angle( (xest'*x) / (x'*x) );
xest = xest*exp(1i*phaseOffset);
errorx = 10*log10(norm(xest - x)^2/norm(x)^2);



[delta SNR test errorx toc]
Errortest(test) = errorx;
runtime(test) = toc;
end
signalnoiseratio(s) = SNR;
Error1(s,nc) = mean(Errortest);
runtime1(s,nc) = mean(runtime);
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


function[f] = innerprod(x,y)
conjy = conj(y);
f = sum(x(:).*conjy(:));
end

function[f] = reversal(x)
f = circshift(flip(x),1);
end