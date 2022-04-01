% %-------------------------------------------------------------------------
% %                       DM FLIGHT MECHANIC
% %-------------------------------------------------------------------------
% %
% %   Full name : Quentin Gourier
% %   Class : A4 SET1
% %   Assignment code : Aé421 (P. Yazigi)
% %
% %-------------------------------------------------------------------------
% %                             DATA
% %-------------------------------------------------------------------------

g = 9.81; %gravity

%aircraft data of lateral stability
h = 6096; %alt
rho = 0.653; %air density
U0 = 94.8; %speed
M = 0.3; %mach number
center_gravity = 0.263;
theta_0 = deg2rad(2.2); % 2.2 deg
omega_0 = 0;

%aircraft data
S = 21.81; %wing area
b = 11.44; %wingspan
c_bar = 2.05; %wing mean chord
A = 6.00; %aspect ratio
e = 0.89; %oswald number
m = 6214; %mass
I_xx = 30048.7; 
I_yy = 26639.8;
I_zz = 55299.8;
I_xz = 606;

C_L_0 = 0.813;
C_D_0 = 0.135;
C_Tx0 = 0.025;
C_m0 = 0;
C_mt0 = 0;

C_m_u = 0;
C_m_alpha = -0.401;
C_m_d_alpha = -5;
C_m_q = -10;
C_m_t_u = 0;
C_m_t_alpha = 0;
C_L_u = 0;
C_L_alpha = 5.22;
C_l_d_alpha = 1.74;
C_L_delta_t = 0;
C_L_q = 3.9;
C_D_u = 0;
C_D_alpha = 0.54;
C_D_d_alpha = 0;
C_D_delta_t = 0;
C_D_q = 0;
C_T_u = 0;
C_L_delta_e = 0.34;
C_D_delta_e = 0;
C_m_delta_e = -0.89;
C_m_delta_t = 0;

C_y_beta = -0.72; C_y_p = 0; C_y_r = 0;
C_l_beta = -0.127; C_l_p = -0.57; C_l_r = 0.2;
C_n_beta = 0.049; C_n_p = -0.045; C_n_r = -0.16;
C_y_delta_r = 0.17; C_l_delta_r = -0.002; C_n_delta_r = -0.073;
C_y_delta_a = 0; C_l_delta_a = 0.14; C_n_delta_a = -0.009;

q_bar = 0.5*rho*U0^2; %dynamic pressure

Yv = (q_bar*S*C_y_beta)/(m*U0);
Yp = (q_bar*S*b*C_y_p)/(2*m*U0);
Yr = (q_bar*S*b*C_y_r)/(2*m*U0);

Lv = (q_bar*S*b*C_l_beta)/(I_xx*U0);
Lp = (q_bar*S*b^2*C_l_p)/(2*I_xx*U0);
Lr = (q_bar*S*b^2*C_l_r)/(2*I_xx*U0);

Nv = (q_bar*S*b*C_n_beta)/(I_zz*U0);
Np = (q_bar*S*b^2*C_n_p)/(2*I_zz*U0);
Nr = (q_bar*S*b^2*C_n_r)/(2*I_zz*U0);

Y_delta_r = (q_bar*S*C_y_delta_r)/(m*U0);
L_delta_r = (q_bar*S*b*C_l_delta_r)/(I_xx*U0);
N_delta_r = (q_bar*S*b*C_n_delta_r)/(I_zz*U0);

Y_delta_a = (q_bar*S*C_y_delta_a)/(m*U0);
L_delta_a = (q_bar*S*b*C_l_delta_a)/(I_xx*U0);
N_delta_a = (q_bar*S*b*C_n_delta_a)/(I_zz*U0);

% aircraft matrix elements of longitudinal stability
q_bar = 0.5*rho*U0^2;
X_u = -((0.5*rho*U0^2*S)/(m*U0))*(2*C_D_0+C_D_u);
X_w = ((q_bar*S)/(m*U0))*(C_L_0-(2/(pi*e*A))*(C_L_0*C_L_alpha));
Z_u = -((0.5*rho*U0^2*S)/(m*U0))*(2*C_L_0+C_L_u);
Z_w = -((0.5*rho*U0^2*S)/(m*U0))*(C_D_0+C_L_alpha);
Z_q = ((0.5*rho*U0^2*S*c_bar)/(2*m*U0))*C_L_q;
M_u = ((0.5*rho*U0^2*S*c_bar)/(I_yy*U0))*C_m_u;
M_w = ((0.5*rho*U0^2*S*c_bar)/(I_yy*U0))*C_m_alpha;
M_d_w = ((0.5*rho*U0^2*S*c_bar^2)/(2*I_yy*U0^2))*C_m_d_alpha;
M_q = ((0.5*rho*U0^2*S*c_bar^2)/(2*I_yy*U0))*C_m_q;

% elevator derivatives
X_delta_e = ((0.5*rho*U0^2*S)/(m*U0))*C_D_delta_e;
Z_delta_e = ((0.5*rho*U0^2*S)/(m*U0))*C_L_delta_e;
M_delta_e = ((0.5*rho*U0^2*S*c_bar)/(I_yy*U0))*C_m_delta_e;

% throttle derivatives
X_delta_t = ((0.5*rho*U0^2*S)/(m*U0))*C_D_delta_t;
Z_delta_t = ((0.5*rho*U0^2*S)/(m*U0))*C_L_delta_t;
M_delta_t = ((0.5*rho*U0^2*S*c_bar)/(I_yy*U0))*C_m_delta_t;



% %-------------------------------------------------------------------------
% %                   LONGITUDINAL STABILITY PART 
% %-------------------------------------------------------------------------

% 1. Equations of longitudinal motion ?
% See the values above and write it on.


% 2. The matrix A of aircraft ?

A = [X_u X_w 0 -g*cos(theta_0);
     Z_u Z_w U0 -g*sin(theta_0);
     M_u+Z_u*M_d_w M_w+Z_w*M_d_w M_q+U0*M_d_w 0;
     0 0 1 0]

B = [X_delta_e X_delta_t;
     Z_delta_e Z_delta_t;
     M_delta_e+Z_delta_e*M_d_w M_delta_t+Z_delta_t*M_d_w;
     0 0]
 
% 3. The characteristic equation ?
% See the matrices and write it on.

% 4. The eigenvalues of the system ?
lambda_vect = eig(A); %vecteur propre
disp(lambda_vect);

lambda_1 = lambda_vect(1,1);
lambda_2 = lambda_vect(2,1);
lambda_3 = lambda_vect(3,1);
lambda_4 = lambda_vect(4,1);

omega = [sqrt(real(lambda_1)^2+imag(lambda_1)^2) sqrt(real(lambda_3)^2+imag(lambda_3)^2)];
omega_n_sp = max(omega)
omega_n_p = min(omega)

ksi = [-real(lambda_1)/omega(1) -real(lambda_3)/omega(2)];
ksi_sp = max(ksi)
ksi_p = min(ksi)

t = linspace(0,400);
t2 = linspace(0,10);

u = exp(-omega_n_p*ksi_p*t).*(U0/10*cos(omega_n_p*sqrt(1-ksi_p^2)*t)+10*omega_n_p*ksi_p/(omega_n_p*sqrt(1-ksi_p^2))*sin(omega_n_p*sqrt(1-ksi_p^2)*t));
theta = exp(-omega_n_p*ksi_p*t).*(theta_0/10*cos(omega_n_p*sqrt(1-ksi_p^2)*t)+theta_0*ksi_p*omega_n_p/(omega_n_p*sqrt(1-ksi_p^2))*sin(omega_n_p*sqrt(1-ksi_p^2)*t));
w = exp(-omega_n_sp*ksi_sp*t2).*(omega_0/10*cos(omega_n_sp*sqrt(1-ksi_sp^2)*t2)+10*ksi_sp*omega_n_sp/(omega_n_sp*sqrt(1-ksi_sp^2))*sin(omega_n_sp*sqrt(1-ksi_sp^2)*t2));
q = exp(-omega_n_sp*ksi_sp*t2).*(theta_0/10*cos(omega_n_sp*sqrt(1-ksi_sp^2)*t2)+theta_0*ksi_sp*omega_n_sp/(omega_n_sp*sqrt(1-ksi_sp^2))*sin(omega_n_sp*sqrt(1-ksi_sp^2)*t2));

% 5. Different modes of longitudinal stability
% a. Short period mode
% b. Phugoid mode


% 6. Curves of longitudinal motion
% a. Axial velocity i function of time
figure
plot(t, u);
legend({'u(t)'})
title('Axial velocity impulse-time response');
xlabel("Time (s)");
ylabel("Axial velocity");
grid("minor");

% b. Angle of attack
figure
plot(t2, w);
legend({'w(t)'})
title('Angle of attack impulse-time response');
xlabel("Time (s)");
ylabel("Angle of attack");
grid("minor");

% c. Pitch rate
figure
plot(t2, q);
legend({'q(t)'})
title('Pitch rate impulse-time response');
xlabel("Time (s)");
ylabel("Pitch rate");
grid("minor");

% d. Pitch angle
figure
plot(t, theta);
legend({'theta(t)'})
title('Pitch angle impulse-time response');
xlabel("Time (s)");
ylabel("Pitch angle");
grid("minor");



% 7. TFs of each variable

s = tf('s');
I = eye(size(A));
TF = zpk((s*I-A)\B)
figure
bode(TF(1,1))
figure
bode(TF(2,1))
figure
bode(TF(3,1))
figure
bode(TF(4,1))




% %-------------------------------------------------------------------------
% %                   LATERAL STABILITY PART 
% %-------------------------------------------------------------------------

% 1. Equations of longitudinal motion ?
% See the values above and write it on.

% 2. The matrix A of aircraft ?

A = [Yv Yp -(U0-Yr) g*cos(theta_0);
     Lv Lp Lr 0;
     Nv Np Nr 0;
     0 1 0 0]

B = [Y_delta_r Y_delta_a;
     L_delta_r L_delta_a;
     N_delta_r N_delta_a
     0 0]

% 3. The characteristic equation ?
% See the matrices and write it on.

% 4. The eigenvalues of the system ?
[vect_propre, val_propre] = eig(A) %vecteur propre

lambda_dutch_roll_1 = val_propre(3,3)
lambda_dutch_roll_2 = val_propre(2,2)
lambda_roll = val_propre(1,1)
lambda_spiral = val_propre(4,4)
 
omega_n_dutch_roll = sqrt(real(lambda_dutch_roll_1)^2+imag(lambda_dutch_roll_1)^2)
omega_n_dutch_roll = sqrt(lambda_dutch_roll_1*lambda_dutch_roll_2)
ksi_dutch_roll = -real(lambda_dutch_roll_1)/omega_n_dutch_roll
ksi_dutch_roll = (lambda_dutch_roll_1+lambda_dutch_roll_2)/(2*omega_n_dutch_roll)

% omega_n_dutch_roll = sqrt(Yv*Nr+Nv*(U0-Yr));
% ksi_dutch_roll = -(Yv+Nr)/(2*omega_n_dutch_roll);

% 5. Different modes of lateral stability
% omega and ksi

% 6. Curves of longitudinal motion
% a. spiral mode
t2 = linspace(0,350);
figure; hold on;
sm1 = vect_propre(1,4)*exp(lambda_spiral*t2);
sm2 = vect_propre(2,4)*exp(lambda_spiral*t2);
sm3 = vect_propre(3,4)*exp(lambda_spiral*t2);
sm4 = vect_propre(4,4)*exp(lambda_spiral*t2);
plot(t2, sm1)
plot(t2, sm2)
plot(t2, sm3)
plot(t2, sm4)
legend({'roll angle', 'yaw rate', 'side slip', 'roll rate'});
title('Spiral mode');
grid("minor");
xlabel("time (seconds)");

% b. rolling mode
t1 = linspace(0,10);
figure; hold on;
% xlim([0 10]);
% ylim([-1 1]);
rm1 = vect_propre(1,1)*exp(lambda_roll*t1);
rm2 = vect_propre(2,1)*exp(lambda_roll*t1);
rm3 = vect_propre(3,1)*exp(lambda_roll*t1);
rm4 = vect_propre(4,1)*exp(lambda_roll*t1);
plot(t1,rm1)
plot(t1,rm2)
plot(t1,rm3)
plot(t1,rm4)
legend({'side velocity', 'roll rate', 'yaw rate', 'side angle'})
title('Rolling mode');
grid("minor");
xlabel("time (seconds)");

% c. dutch roll mode
t3 = linspace(0,50);
figure; hold on;
% xlim([0 50]);
% ylim([-1 1]);
dr1 = vect_propre(1,2)*(exp(real(lambda_dutch_roll_2).*t3).*(cos(omega_n_dutch_roll.*t3)+ksi_dutch_roll*sin(omega_n_dutch_roll.*t3)));
dr2 = vect_propre(2,2)*(exp(real(lambda_dutch_roll_2).*t3).*(cos(omega_n_dutch_roll.*t3)+ksi_dutch_roll*sin(omega_n_dutch_roll.*t3)));
dr3 = vect_propre(3,2)*(exp(real(lambda_dutch_roll_2).*t3).*(cos(omega_n_dutch_roll.*t3)+ksi_dutch_roll*sin(omega_n_dutch_roll.*t3)));
dr4 = vect_propre(4,2)*(exp(real(lambda_dutch_roll_2).*t3).*(cos(omega_n_dutch_roll.*t3)+ksi_dutch_roll*sin(omega_n_dutch_roll.*t3)));
plot(t3, dr1)
plot(t3, dr2)
plot(t3, dr3)
plot(t3, dr4)
legend({'side velocity', 'roll rate', 'yaw rate', 'side angle'});
title('Dutch roll mode');
grid("minor");
xlabel("time (seconds)");

% Stability balance data part
% We precise these values have been choosen with adjustment
new_omega_n_dutch_roll = 1.2;
new_ksi_dutch_roll = 0.24;
k1 = -284.6;
k2 = 6.45;
new_lambda_dutch_roll_1 = -0.288 -1.165i;
new_lambda_dutch_roll_2 = -0.288 +1.165i;

t3= linspace(0,20);
figure; hold on;
new_dr1 = vect_propre(1,2).*exp(real(new_lambda_dutch_roll_1)*t3).*(cos(new_omega_n_dutch_roll.*t3)+new_ksi_dutch_roll*sin(new_omega_n_dutch_roll.*t3));
new_dr2 = vect_propre(2,2).*exp(real(new_lambda_dutch_roll_1)*t3).*(cos(new_omega_n_dutch_roll.*t3)+new_ksi_dutch_roll*sin(new_omega_n_dutch_roll.*t3));
new_dr3 = vect_propre(3,2).*exp(real(new_lambda_dutch_roll_1)*t3).*(cos(new_omega_n_dutch_roll.*t3)+new_ksi_dutch_roll*sin(new_omega_n_dutch_roll.*t3));
new_dr4 = vect_propre(4,2).*exp(real(new_lambda_dutch_roll_1)*t3).*(cos(new_omega_n_dutch_roll.*t3)+new_ksi_dutch_roll*sin(new_omega_n_dutch_roll.*t3));
plot(t3, new_dr1)
plot(t3, new_dr2)
plot(t3, new_dr3)
plot(t3, new_dr4)
legend({'side velocity', 'roll rate', 'yaw rate', 'side angle'});
title('Balanced dutch roll mode');
grid("minor");
xlabel("time (seconds)");

% 7. TFs of each variable
s = tf('s');
I = eye(size(A));
TF = zpk((s*I-A)\B);
figure
bode(TF(1,1))
figure
bode(TF(2,1))
figure
bode(TF(3,1))
figure
bode(TF(4,1))
figure
bode(TF(1,2))
figure
bode(TF(2,2))
figure
bode(TF(3,2))
figure
bode(TF(4,2))


%END OF TP
