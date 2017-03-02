%% unscented trafo test

% posterior state

T = 1;
x = zeros(4,1);

F = [1 0 T 0
     0 1 0 T
     0 0 1 0
     0 0 0 1];
 
P = eye(4);

var_v1 = 0.1;
var_v2 = 0.1;
var_w1 = 0.2;
var_w2 = 0.2;

G = [T^2/2 0 
     0     T^2/2
     T     0
     0     T];
Q = G* [var_v1 0 ; 0 var_v2 ]* G'

% linear prediction
xpred = F*x
Ppred = F*P*F' + Q

% unscented prediction state and process noise
x_aug = [x ; 0 ; 0];
P_aug = [P zeros(4,2) ; zeros(2,4) [var_v1 0 ; 0 var_v2 ]];

n_aug = 6;
kappa = 3-n_aug;

w = zeros(2*n_aug+1);
w(1) = kappa/(kappa+n_aug);
for i=2:13
  w(i) = 0.5/(n_aug+kappa);
end

Xi_aug = zeros(n_aug,2*n_aug+1);
sP_aug = chol(P_aug);

Xi_aug(:,1) = x_aug;
Xi_aug(:,2) = x_aug + sqrt(n_aug+kappa) * sP_aug(:,1);
Xi_aug(:,3) = x_aug - sqrt(n_aug+kappa) * sP_aug(:,1);
Xi_aug(:,4) = x_aug + sqrt(n_aug+kappa) * sP_aug(:,2);
Xi_aug(:,5) = x_aug - sqrt(n_aug+kappa) * sP_aug(:,2);
Xi_aug(:,6) = x_aug + sqrt(n_aug+kappa) * sP_aug(:,3);
Xi_aug(:,7) = x_aug - sqrt(n_aug+kappa) * sP_aug(:,3);
Xi_aug(:,8) = x_aug + sqrt(n_aug+kappa) * sP_aug(:,4);
Xi_aug(:,9) = x_aug - sqrt(n_aug+kappa) * sP_aug(:,4);
Xi_aug(:,10) = x_aug + sqrt(n_aug+kappa) * sP_aug(:,5);
Xi_aug(:,11) = x_aug - sqrt(n_aug+kappa) * sP_aug(:,5);
Xi_aug(:,12) = x_aug + sqrt(n_aug+kappa) * sP_aug(:,6);
Xi_aug(:,13) = x_aug - sqrt(n_aug+kappa) * sP_aug(:,6);

xpred_un = 0;
Ppred_un = zeros(4,4);

% prediction
Xi_pred = F * Xi_aug(1:4,:) + G*Xi_aug(5:6,:);

for i=1:2*n_aug+1
    xpred_un = xpred_un + w(i)* Xi_pred(:,i);
end

for i=1:2*n_aug+1
    Ppred_un = Ppred_un + w(i)* (Xi_pred(:,i) - xpred_un)*(Xi_pred(:,i) - xpred_un)';
end

xpred_un
Ppred_un
