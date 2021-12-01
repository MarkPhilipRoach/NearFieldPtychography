clear all
Tests = 100;
ca = cell(1,4);
%% Assigning variables

d = 315;
Xtilde = zeros(d,d);
Weight = zeros(d,d);
Dmatrix = zeros(d,d);
objectX = randn(d,Tests)+ 1i*randn(d,Tests);

signalnoiseratio = zeros(1,4);
Error1 = zeros(4,4);
runtime1 = zeros(1,1);
nc = 0;
for delta = [2 5 8 11]
nc = nc + 1; 
ca{nc} = sprintf('delta = %d', delta);
D = d*(2*delta-1);
shifts = d;
K = 0:shifts-1; %set of shifts
L = 0:2*delta-2;% %set of frequencies
Knum = length(K);
Lnum = length(L);
Noise = randn(Tests*Knum,Lnum) + 1i*randn(Tests*Knum,Lnum);

a = max(4,(delta-1)/2);
pointspread = zeros(d,1);
maskorg = zeros(d,1);
for n = 1:d
    pointspread(n) = exp(-2*pi*1i*n^2/(2*delta-1));
end
for n = 1:delta
    maskorg(n) = (exp((-n+1)/a))/((2*delta-1)^(1/4))*exp(2*pi*1i*n^2/(2*delta-1));
end




%% Constructing block measurment matrix M'
s = 0;
for SNR = 20:20:80
s = s + 1;
Errortest = zeros(Tests,1);
runtime = zeros(Tests,1);
for test = 1:Tests
object = objectX(:,test);
tic



matrixmask = zeros(Lnum,d);
for l = 1:Lnum
   shiftflipp = circshift(reversal(pointspread),-l+1);
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
maskmatrix = repmat(masklmatrix,Knum,1);
for i=1:Knum
    maskmatrix((i-1)*(2*delta-1)+1:i*(2*delta-1),:) = circshift(maskmatrix((i-1)*(2*delta-1)+1:i*(2*delta-1),:) ,(i-1)*(2*delta-1),2);
end



%% Constructing measurements
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
Yb = zeros(Knum,Lnum);
for k = 1:Knum
    for l = 1:Lnum
         Yb(k,l) = abs(innerprod(matrixmask(l,:),circshift(object,-k+1)))^2;
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

circy = circshift(y,delta-1);
X = zeros(d,d);
for n=1:d
    X(:,n) = circshift([circy((n*(2*delta-1)-(2*delta-2)):n*(2*delta-1)); zeros(d-(2*delta-1),1)],n-delta);
end
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
phaseOffset = angle( (xest'*object) / (object'*object) );
xest = xest*exp(1i*phaseOffset);
errorx = 10*log10(norm(xest - object)^2/norm(object)^2);
[delta SNR test errorx toc]
Errortest(test) = errorx;
runtime(test) = toc;
end
signalnoiseratio(s) = SNR;
Error1(s,nc) = mean(Errortest);
runtime1(s,nc) = mean(runtime);
end
end
plot(signalnoiseratio,Error1(:,1),'-b','LineWidth',1.5)
hold on
plot(signalnoiseratio,Error1(:,2),'--r','LineWidth',1.5)
hold on
plot(signalnoiseratio,Error1(:,3),':g','LineWidth',1.5)
hold on
plot(signalnoiseratio,Error1(:,4),'-.k','LineWidth',1.5)

xlabel({'SNR (in dB)'})
ylabel({'Reconstruction Error (in dB)'})
title({'SNR vs Reconstruction Error'})
xticks(20:10:80) 
legend(ca, 'Location', 'northeast')


figure()

X = categorical({'2','5','8','11'});
X = reordercats(X,{'2','5','8','11'});
Y = [mean(runtime1(:,1)) mean(runtime1(:,2)) mean(runtime1(:,3)) mean(runtime1(:,4))];
plot(X,Y,'Marker','o','Color', 'b','LineWidth',2)

xlabel({'\delta'})
ylabel({'Runtime (in seconds)'})
title({'\delta vs Runtime'})

%% Pre-Assigned Functions
function[f] = reversal(x)
f = circshift(flip(x),1);
end

function[f] = innerprod(x,y)
conjy = conj(y);
f = sum(x(:).*conjy(:));
end
