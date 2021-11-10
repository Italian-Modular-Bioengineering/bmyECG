clear classes;

function [K X Y P] = EKFstep(ES, R, Q, Yd, X, P)
Ax = ES.Ad(X);
P = Ax*P*Ax' + Q;
X = ES.fd(X);
Cx = ES.Cd(X);
K = P*Cx'*(Cx*P*Cx' + R)^-1;
P = P - K*Cx*P;
Y = ES.g(X);
corr = K*(Yd-Y);
X = X + corr;
X(4) = max(X(4), ES.minPulse); %frequency saturation

endfunction


SR = 2000;
dt = 1/SR;

%ES = ECG_SYS(dt);
ES = ECG_SYS_euler(dt);

noise = 1e-2;

T = 10;
L = T*SR;
t = (0:L-1)*dt;
ym = zeros(1,L);
omegam = 2*pi*2;
Xm0 = [0 1 0 omegam]';
Xm = zeros(4,L);
Xm(:,1) = Xm0;

# imposing initial values
omega0 = 2*pi*1;
X0 = [ 0.1 0.2 0 omega0 ]';
#initializing null vector which represent the state
X = zeros(4,L);
# assegnation of the initial value for the state
X(:,1) = X0;
#output vector initialization
y = zeros(1,L);

# initializing matrices
R = 10;
Q = diag([ 0.1 0.1 0 1000 ]);
P = Q;

for i = 1:L-1
  ym(i) = ES.g(Xm(:,i)) + randn(1,1)*noise;
  Xm(:,i+1) = ES.fd(Xm(:,i));
  
  [ K X(:,i+1) y(i) P] = EKFstep(ES, R, Q, ym(i), X(:,i), P);
  
end

figure(1); plot(t, ym, t, y); legend;
figure(2); plot(X'); legend;
figure(3); plot(t, X(4,:)');
  