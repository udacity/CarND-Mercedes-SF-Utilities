%% unscented trafo test

% posterior state
x = 1;
F = 1.2;
P = 0.04;
Q = 0.1;

% linear prediction
xpred = F*x
Ppred = F*P*F' + Q

% unscented prediction state uncertainty only
n = 1; 
kappa = 3-n;

w0 = kappa/(kappa+n);
w1 = 0.5/(n+kappa);
w2 = 0.5/(n+kappa);

% w0 = 1/3;
% w1 = 1/3;
% w2 = 1/3;

xi0 = x;
xi1 = x + sqrt(n+kappa) * sqrt(P);
xi2 = x - sqrt(n+kappa) * sqrt(P);

xi0pred = F*xi0;
xi1pred = F*xi1;
xi2pred = F*xi2;

xpred_u = w0*xi0pred + w1*xi1pred + w2*xi2pred
Ppred_u = w0*(xi0pred-xpred_u)^2 + w1*(xi1pred-xpred_u)^2 +w2*(xi2pred-xpred_u)^2 + Q


% unscented prediction state and process noise
x_aug = [1 ; 0];
P_aug = [P 0 ; 0 Q];

n = 2;
kappa = 3-n;

w = zeros(2*n+1);
w(1) = kappa/(kappa+n);
w(2) = 0.5/(n+kappa);
w(3) = 0.5/(n+kappa);
w(4) = 0.5/(n+kappa);
w(5) = 0.5/(n+kappa);

Xi_aug = zeros(n,2*n+1);
sP_aug = chol(P_aug);

Xi_aug(:,1) = x_aug;
Xi_aug(:,2) = x_aug + sqrt(n+kappa) * sP_aug(:,1);
Xi_aug(:,3) = x_aug - sqrt(n+kappa) * sP_aug(:,1);
Xi_aug(:,4) = x_aug + sqrt(n+kappa) * sP_aug(:,2);
Xi_aug(:,5) = x_aug - sqrt(n+kappa) * sP_aug(:,2);


xpred_un = 0;
Ppred_un = 0;


% prediction
Xi_pred = F * Xi_aug(1,:) + Xi_aug(2,:);

for i=1:5
    xpred_un = xpred_un + w(i)* Xi_pred(1,i);
end

for i=1:5
    Ppred_un = Ppred_un + w(i)* (Xi_pred(1,i) - xpred_un)^2;
end

xpred_un
Ppred_un
