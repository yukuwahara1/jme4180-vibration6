close all
clc
E = 193e9;       % young's modulus (Pa)
rho = 8000;     % density (kg/m^3)
L = .599;        % length (m)
b = 0.019;      % width (m)
h = 0.0063;     % thickness (m)
I  =  b*h^3/12; % area moment of inertia
A = b*h;        % cross section area
B1L = 1.875;    % eignevalue 1
B2L = 4.694;    % eigenvalue 2
B3L = 7.855;    % eigenvalue 3

% COMPUTE NATURAL FREQUENCIES
BL = [B1L;B2L;B3L];                    
wn = (BL.^2)*sqrt(E*I/(rho*A*L^4));     % natural frequencies in rad/sec
fn = wn/(2*pi)     % natural frequencies in Hz

% COMPUTE MODE SHAPES
alpha = (sin(BL) + sinh(BL))./(cos(BL) + cosh(BL));
dx = L/20;
x=0:dx:L;
W1 =  (sin(BL(1)*x/L) - sinh(BL(1)*x/L) -alpha(1)*(cos(BL(1)*x/L) - cosh(BL(1)*x/L)));
W2 =  (sin(BL(2)*x/L) - sinh(BL(2)*x/L) -alpha(2)*(cos(BL(2)*x/L) - cosh(BL(2)*x/L)));
W3 =  (sin(BL(3)*x/L) - sinh(BL(3)*x/L) -alpha(3)*(cos(BL(3)*x/L) - cosh(BL(3)*x/L)));
% SCALE SO W(L)=1
W1  = W1/W1(end);
W2  = W2/W2(end);
W3  = W3/W3(end);

% MASS NORMALIZATION OF MODE SHAPES
m1 = sum(rho*A*W1.^2*dx);
W1  = W1/sqrt(m1);
mx1 = 1.2*max(abs(W1));
m2 = sum(rho*A*W2.^2*dx);
W2  = W2/sqrt(m2);
mx3 = 1.2*max(abs(W3));
m3 = sum(rho*A*W3.^2*dx);
W3  = W3/sqrt(m3);
mx2 = 1.2*max(abs(W2));
Ba1L = .07198;    % eignevalue 1
Ba2L = 1.894;    % eigenvalue 2
Ba3L = 3.026;    % eigenvalue 3

% COMPUTE NATURAL FREQUENCIES
BcL = [Ba1L;Ba2L;Ba3L];                    
wnb = (BcL.^2)*sqrt(E*I/(rho*A*L^4));     % natural frequencies in rad/sec
fnb = wnb/(2*pi)                         % natural frequencies in Hz

% COMPUTE MODE SHAPES
alphab = (sin(BcL) + sinh(BcL))./(cos(BcL) + cosh(BcL));
dxb = L/20;
xb=0:dxb:L;
W1a =  (sin(BcL(1)*xb/L) - sinh(BcL(1)*xb/L) -alphab(1)*(cos(BcL(1)*xb/L) - cosh(BcL(1)*xb/L)));
W2a =  (sin(BcL(2)*xb/L) - sinh(BcL(2)*xb/L) -alphab(2)*(cos(BcL(2)*xb/L) - cosh(BcL(2)*xb/L)));
W3a =  (sin(BcL(3)*xb/L) - sinh(BcL(3)*xb/L) -alphab(3)*(cos(BcL(3)*xb/L) - cosh(BcL(3)*xb/L)));
% SCALE SO W(L)=1
W1a  = W1a/W1a(end);
W2a  = W2a/W2a(end);
W3a  = W3a/W3a(end);
% MASS NORMALIZATION OF MODE SHAPES
m1a = sum(rho*A*W1a.^2*dxb);
W1a  = W1a/sqrt(m1a);
mx1a = 1.2*max(abs(W1a));
m2a = sum(rho*A*W2a.^2*dxb);
W2a  = W2a/sqrt(m2a);
mx3a = 1.2*max(abs(W3a));
m3a = sum(rho*A*W3a.^2*dxb);
W3a  = W3a/sqrt(m3a);
mx2a = 1.2*max(abs(W2a));
% solidworks simulation data
Bb1L = 0.11985;    % eignevalue 1, need to be changed
Bb2L = 0.7507;    % eigenvalue 2, need to be changed
Bb3L = 2.1003;    % eigenvalue 3, need to be changed
% sw = SolidWorks
% COMPUTE NATURAL FREQUENCIES
BswL = [Bsw1L;Bsw2L;Bsw3L];                    
wnsw = (BswL.^2)*sqrt(E*I/(rho*A*L^4));     % natural frequencies in rad/sec
fnsw = wnsw/(2*pi)                         % natural frequencies in Hz
% COMPUTE MODE SHAPES
alphab = (sin(BswL) + sinh(BswL))./(cos(BswL) + cosh(BswL));
dxb = L/20;
xb=0:dxb:L;
W1sw =  (sin(BswL(1)*xb/L) - sinh(BswL(1)*xb/L) -alphab(1)*(cos(BswL(1)*xb/L) - cosh(BswL(1)*xb/L)));
W2sw =  (sin(BswL(2)*xb/L) - sinh(BswL(2)*xb/L) -alphab(2)*(cos(BswL(2)*xb/L) - cosh(BswL(2)*xb/L)));
W3sw =  (sin(BswL(3)*xb/L) - sinh(BswL(3)*xb/L) -alphab(3)*(cos(BswL(3)*xb/L) - cosh(BswL(3)*xb/L)));
% SCALE SO W(L)=1
W1sw = W1sw/W1sw(end);
W2sw  = W2sw/W2sw(end);
W3sw  = W3sw/W3sw(end);
% MASS NORMALIZATION OF MODE SHAPES
m1sw = sum(rho*A*W1sw.^2*dxb);
W1sw  = W1sw/sqrt(m1sw);
mx1sw = 1.2*max(abs(W1sw));
m2sw = sum(rho*A*W2sw.^2*dxb);
W2sw  = W2sw/sqrt(m2sw);
mx3sw = 1.2*max(abs(W3sw));
m3sw = sum(rho*A*W3sw.^2*dxb);
W3sw  = W3sw/sqrt(m3sw);
mx2sw = 1.2*max(abs(W2sw));
% PLOT NORMALIZED MODE SHAPES FROM MATLAB
% Note the different length axis options! 
figure 
plot(100*x,W1); 
title('Mode 1')
xlabel(sprintf('x (cm) \n'))
ylabel('Amplitude (Normalized)')
hold on 
plot (100*xb,W1a);
hold on 
plot
hold off
figure 
plot(x,W2); 
title('Mode 2')
xlabel(sprintf('x (m) \n'))
ylabel('Amplitude (Normalized)')
hold on
plot(xb, W2a);
hold on
plot
hold off
figure 
plot(x/L(end),W3);
title('Mode 3')
xlabel('x (Normalized)')
ylabel('Amplitude (Normalized)')
hold on
plot(xb/L(end),W3a);
hold on
plot(
hold off
