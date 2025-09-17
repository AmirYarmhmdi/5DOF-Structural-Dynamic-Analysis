%% Amir Yarmohamadi 9411182
%   MATLAB script for 5-DOF shear building modal analysis
% - Builds M and K
% - Eigen-solution, modal properties, Rayleigh damping
% - Step-force and harmonic-force modal time histories
% - Base shear / overturning moment / drifts
% - El-Centro response history (modal superposition, central difference–style)
% - Response spectrum–style envelope (using max modal responses)
clc; clear; close all;

%% ---------------------------
%% Parameters
K    = 315;         % Story stiffness [kips/in]
M    = 135;         % Lumped floor mass [kips]
zeta = 0.03;        % Target modal damping ratio for Rayleigh construction (approx.)
h    = [12 24 36 48 60]; % Floor elevations [ft] (for resultants like overturning)
H    = 1:5;         % Story numbers for plotting (y-axis)

nDOF = 5;

%% ---------------------------
%% Assemble stiffness and mass matrices
k = zeros(nDOF);    % 5x5 stiffness
m = zeros(nDOF);    % 5x5 mass

% Tridiagonal shear-building stiffness (fixed base)
for o = 1:(nDOF-1)
    k(o,  o)   = k(o,  o)   + 2*K;
    k(o,  o+1) = k(o,  o+1) - K;
    k(o+1,o)   = k(o+1,o)   - K;
    k(o+1,o+1) = k(o+1,o+1) + K;
end

% Lumped mass on each floor
for i = 1:nDOF
    m(i,i) = M;
end

%% ---------------------------
%% Eigen-solution (generalized)
[phi, w] = eig(k, m);   % k*phi = w*m*phi
W  = diag(w);           % squared circular frequencies
Wn = sqrt(W);           % circular frequencies [rad/s]

% Mode shape columns
phi1 = phi(:,1); phi2 = phi(:,2); phi3 = phi(:,3); phi4 = phi(:,4); phi5 = phi(:,5);

% Frequencies and periods
Frequency = zeros(1,nDOF);
Period    = zeros(1,nDOF);
for i = 1:nDOF
    Frequency(i) = Wn(i)/(2*pi); % [Hz]
    Period(i)    = 1/Frequency(i);
end

%% ---------------------------
%% Rayleigh damping matrix (using modes 1 and 3 as targets)
a0 = 2*zeta * ((2*Wn(1)*Wn(3))/(Wn(1)+Wn(3))); % mass-proportional coeff
a1 = 2*zeta * (1/(Wn(1)+Wn(3)));               % stiffness-proportional coeff

zeta_i = zeros(1,nDOF); % effective modal damping ratios from Rayleigh C
for i = 1:nDOF
    zeta_i(i) = (a0/(2*Wn(i))) + (a1*Wn(i))/2;
end

c = a0.*m + a1.*k; % Rayleigh damping matrix

%% ---------------------------
%% Display basics
disp('***** Answers of Problems 1 & 2 *****');
disp('M =');  disp(m);
disp('K =');  disp(k);
disp('C =');  disp(c);
disp('Wn (rad/s) ='); disp(Wn');
disp('Modal damping ratios from Rayleigh (zeta_i) ='); disp(zeta_i);
disp('Mode shapes (phi) ='); disp(phi);

%% ---------------------------
%% Plot mode shapes vs story index
figure(1);
for i = 1:nDOF
    subplot(1,nDOF,i);
    plot(phi(:,i), H, '-o'); grid on;
    xlabel(['\phi_',num2str(i)]); ylabel('Story');
    title(['Mode ',num2str(i)]);
end

%% ================================================================
%% MODAL ANALYSIS for step force
%% ================================================================
% Step load vector [kips] (example pattern)
pt = [1000 0 0 500 -500]'; % step applied at t>=0

% Generalized modal quantities
MM = zeros(1,nDOF); KK = zeros(1,nDOF); Pt = zeros(1,nDOF); Wd = zeros(1,nDOF);
for i=1:nDOF
    MM(i) = phi(:,i)'*m*phi(:,i);    % generalized mass
    KK(i) = phi(:,i)'*k*phi(:,i);    % generalized stiffness
    Pt(i) = phi(:,i)'*pt;            % generalized force (step magnitude)
    Wd(i) = Wn(i)*sqrt(1 - zeta_i(i)^2); % damped circular frequency
end

% Time vector
t = 1:0.05:200; % 3981 samples

% Modal coordinates for unit step response under modal force Pt(i)
q = zeros(nDOF, numel(t));
for i=1:nDOF
    q(i,:) = (Pt(i)/KK(i)) * ( 1 ...
              - exp(-zeta_i(i)*Wn(i)*t) .* ( cos(Wd(i)*t) ...
              + (zeta_i(i)/sqrt(1 - zeta_i(i)^2)) .* sin(Wd(i)*t) ) );
end

% Physical displacements (superposition)
Ut = phi(:,1)*q(1,:) + phi(:,2)*q(2,:) + phi(:,3)*q(3,:) + phi(:,4)*q(4,:) + phi(:,5)*q(5,:);

% Story drifts time-history (absolute differences of floors)
Drift1 = Ut(1,:);
Drift2 = Ut(2,:) - Ut(1,:);
Drift3 = Ut(3,:) - Ut(2,:);
Drift4 = Ut(4,:) - Ut(3,:);
Drift5 = Ut(5,:) - Ut(4,:);

driftt    = [Drift1;Drift2;Drift3;Drift4;Drift5];
driftmax  = [max(abs(Drift1));max(abs(Drift2));max(abs(Drift3));max(abs(Drift4));max(abs(Drift5))];

% Base shear and overturning moment (K times drift approximates story shear in shear model)
vb  = zeros(nDOF, numel(t));
for i=1:nDOF
    vb(i,:) = K .* driftt(i,:);
end
Vb = sum(vb, 1);           % base shear time history (sum of story shears)
mb = zeros(nDOF, numel(t));
for i=1:nDOF
    mb(i,:) = h(i) .* vb(i,:);
end
Mb = sum(mb, 1);           % overturning moment time history

% Max responses
maxU = zeros(1,nDOF);
for i=1:nDOF
    maxU(i) = max(abs(Ut(i,:)));
end
V  = driftmax(:) * K;      % peak story shear from peak drift
Vtotal_step = sum(V);      % base shear envelope (sum of peak story shears)
mm = zeros(1,nDOF);
for i=1:nDOF
    mm(i) = V(i) * h(i);
end
Mtotal_step = sum(mm);     % overturning moment envelope

disp('****** Answers of Problem 3 (Step Force) ******');
disp('Roof displacement time-history Ut(5,:):'); disp(Ut(5,:));
disp('Roof drift time-history Drift5:');         disp(Drift5);
disp('Base shear time-history Vb:');              disp(Vb);
disp('Overturning moment time-history Mb:');      disp(Mb);
disp('Base shear envelope (sum of peak story shears):'); disp(Vtotal_step);
disp('OTM envelope (sum(V_i*h_i)) :');                        disp(Mtotal_step);

figure(2)
subplot(2,3,1); plot(t, Ut(5,:));    grid on; xlabel('t [s]'); ylabel('u_{roof} [in]'); title('Roof Displacement');
subplot(2,3,2); plot(t, Drift5);     grid on; xlabel('t [s]'); ylabel('\Delta_{roof} [in]'); title('Roof Drift');
subplot(2,3,3); plot(t, Vb);         grid on; xlabel('t [s]'); ylabel('V_b [kips]'); title('Base Shear (time history)');
subplot(2,3,4); plot(t, Mb);         grid on; xlabel('t [s]'); ylabel('M_b [kip-ft]'); title('OTM (time history)');
subplot(2,3,5); plot(abs(maxU), H);  grid on; xlabel('|u_{max}| [in]'); ylabel('Story'); title('Max Story Displacements');
subplot(2,3,6); plot(abs(driftmax),H); grid on; xlabel('|drift_{max}| [in]'); ylabel('Story'); title('Max Story Drifts');

%% ================================================================
%% MODAL ANALYSIS for harmonic force (example P0 and excitation)
%% ================================================================
p0t = [1000 0 0 500 -500]; % harmonic force amplitude vector [kips]
Ptt = zeros(1,nDOF);
for i=1:nDOF
    Ptt(i) = phi(:,i)' * p0t';
end

% NOTE: The original code used a closed-form transient for harmonic start.
% Keeping the same formula the user provided.
qq = zeros(nDOF, numel(t));
for i=1:nDOF
    qq(i,:) = (Ptt(i)/KK(i)) * (1/(2*zeta_i(i))) * ...
              ( exp(-zeta_i(i)*Wn(i)*t) .* ( cos(Wd(i)*t) + (zeta_i(i)/sqrt(1-zeta_i(i)^2))*sin(Wd(i)*t) ) ...
                - cos(Wn(i)*t) );
end

Utt = phi(:,1)*qq(1,:) + phi(:,2)*qq(2,:) + phi(:,3)*qq(3,:) + phi(:,4)*qq(4,:) + phi(:,5)*qq(5,:);

% Drifts for harmonic case
Drifft1 = Utt(1,:);
Drifft2 = Utt(2,:) - Utt(1,:);
Drifft3 = Utt(3,:) - Utt(2,:);
Drifft4 = Utt(4,:) - Utt(3,:);
Drifft5 = Utt(5,:) - Utt(4,:);

drifttt   = [Drifft1;Drifft2;Drifft3;Drifft4;Drifft5];
drifttmax = [max(abs(Drifft1));max(abs(Drifft2));max(abs(Drifft3));max(abs(Drifft4));max(abs(Drifft5))];

vbb = zeros(nDOF, numel(t));
for i=1:nDOF
    vbb(i,:) = K .* drifttt(i,:);
end
Vbb = sum(vbb,1);
Vtotal_env_harmonic = max(Vbb);  %#ok<NASGU> % peak base shear over time (not used later)

mbb = zeros(nDOF, numel(t));
for i=1:nDOF
    mbb(i,:) = h(i) .* vbb(i,:);
end
Mbb = sum(mbb,1);

maxxU = zeros(1,nDOF);
for i=1:nDOF
    maxxU(i) = max(abs(Utt(i,:)));
end

Vv      = drifttmax(:) * K;
Vtotal  = sum(Vv);       % envelope base shear from peak story drifts
mmm     = zeros(1,nDOF);
for i=1:nDOF
    mmm(i) = Vv(i) * h(i);
end
Mmm = sum(mmm);

disp('****** Answers of Problem 4 (Harmonic) ******');
disp('Roof displacement time-history Utt(5,:):'); disp(Utt(5,:));
disp('Roof drift time-history Drifft5:');         disp(Drifft5);
disp('Base shear time-history Vbb:');              disp(Vbb);
disp('OTM time-history Mbb:');                    disp(Mbb);
disp('Base shear envelope (sum of peak story shears):'); disp(Vtotal);
disp('OTM envelope (sum(V_i*h_i)) :');                        disp(Mmm);

figure(3)
subplot(2,3,1); plot(t, Utt(5,:));    grid on; xlabel('t [s]'); ylabel('u_{roof} [in]'); title('Roof Displacement (harmonic)');
subplot(2,3,2); plot(t, Drifft5);     grid on; xlabel('t [s]'); ylabel('\Delta_{roof} [in]'); title('Roof Drift (harmonic)');
subplot(2,3,3); plot(t, Vbb);         grid on; xlabel('t [s]'); ylabel('V_b [kips]'); title('Base Shear (harmonic)');
subplot(2,3,4); plot(t, Mbb);         grid on; xlabel('t [s]'); ylabel('M_b [kip-ft]'); title('OTM (harmonic)');
subplot(2,3,5); plot(abs(maxxU), H);  grid on; xlabel('|u_{max}| [in]'); ylabel('Story'); title('Max Story Displacements (harmonic)');
subplot(2,3,6); plot(abs(drifttmax),H); grid on; xlabel('|drift_{max}| [in]'); ylabel('Story'); title('Max Story Drifts (harmonic)');

%% ================================================================
%% Response History Analysis for El-Centro (PGA scaling 0.65g)
%% ================================================================
% Modal participation factors and effective heights
Mn = zeros(1,nDOF);
Mn(1) = M*(phi1(1)^2 + phi1(2)^2 + phi1(3)^2 + phi1(4)^2 + phi1(5)^2);
Mn(2) = M*(phi2(1)^2 + phi2(2)^2 + phi2(3)^2 + phi2(4)^2 + phi2(5)^2);
Mn(3) = M*(phi3(1)^2 + phi3(2)^2 + phi3(3)^2 + phi3(4)^2 + phi3(5)^2);
Mn(4) = M*(phi4(1)^2 + phi4(2)^2 + phi4(3)^2 + phi4(4)^2 + phi4(5)^2);
Mn(5) = M*(phi5(1)^2 + phi5(2)^2 + phi5(3)^2 + phi5(4)^2 + phi5(5)^2);

Lh = zeros(1,nDOF);
Lh(1) = M*(phi1(1)+phi1(2)+phi1(3)+phi1(4)+phi1(5));
Lh(2) = M*(phi2(1)+phi2(2)+phi2(3)+phi2(4)+phi2(5));
Lh(3) = M*(phi3(1)+phi3(2)+phi3(3)+phi3(4)+phi3(5));
Lh(4) = M*(phi4(1)+phi4(2)+phi4(3)+phi4(4)+phi4(5));
Lh(5) = M*(phi5(1)+phi5(2)+phi5(3)+phi5(4)+phi5(5));

Gama = zeros(1,nDOF); % participation factors
for i=1:nDOF
    Gama(i) = Lh(i)/Mn(i);
end
disp('Participation factors (Gamma):'); disp(Gama);

S = zeros(nDOF,nDOF); % modal influence vectors (mass*mode scaled by Gamma)
for i=1:nDOF
    S(:,i) = Gama(i) .* m * phi(:,i);
end

Lteta = zeros(1,nDOF); % moment arms for overturning
Lteta(1) = M*(h(1)*phi1(1)+h(2)*phi1(2)+h(3)*phi1(3)+h(4)*phi1(4)+h(5)*phi1(5));
Lteta(2) = M*(h(1)*phi2(1)+h(2)*phi2(2)+h(3)*phi2(3)+h(4)*phi2(4)+h(5)*phi2(5));
Lteta(3) = M*(h(1)*phi3(1)+h(2)*phi3(2)+h(3)*phi3(3)+h(4)*phi3(4)+h(5)*phi3(5));
Lteta(4) = M*(h(1)*phi4(1)+h(2)*phi4(2)+h(3)*phi4(3)+h(4)*phi4(4)+h(5)*phi4(5));
Lteta(5) = M*(h(1)*phi5(1)+h(2)*phi5(2)+h(3)*phi5(3)+h(4)*phi5(4)+h(5)*phi5(5));

Mstr = zeros(1,nDOF); % modal base shear factors
Hstr = zeros(1,nDOF); % effective modal heights
for i=1:nDOF
    Mstr(i) = Gama(i)*Lh(i);
    Hstr(i) = Lteta(i)/Lh(i);
end

%% ---------------------------
%% Response History (modal superposition via central-difference-like recursion)
tt       = 0:0.01:4;    % 401 samples
delta_t  = 0.01;
scalePGA = 0.65;        % scale factor

% Load ground acceleration record (column vector, units consistent with forcing model)
% Put the file next to this script, or adjust path accordingly.
% 'data' should be length >= 401 for this section.
data = load('elcentro.txt'); % <-- place your acceleration record here

% Mode 1
C1     = 2*Wn(1)*MM(1)*zeta_i(1);
U1     = zeros(1, numel(tt)); U1(1) = 0;
K_hat1 = (MM(1)/(delta_t^2)) + (C1/(2*delta_t));
a11    = (MM(1)/delta_t^2) - (C1/(2*delta_t));
b1     = KK(1) - (2*MM(1)/(delta_t^2));
p      = data(:)' * MM(1);    % modal force proportional to ground acc. * modal mass
p_hat1 = zeros(1, numel(tt)); p_hat1(1) = p(1);
U1(2)  = (p_hat1(1)/K_hat1)*scalePGA;
for i=2:400
    p_hat1(i) = p(i) - a11*U1(i-1) - b1*U1(i);
    U1(i+1)   = (p_hat1(i)/K_hat1)*scalePGA;
end

% Mode 2
C2     = 2*Wn(2)*MM(2)*zeta_i(2);
U2     = zeros(1, numel(tt)); U2(1) = 0;
K_hat2 = (MM(2)/(delta_t^2)) + (C2/(2*delta_t));
a2     = (MM(2)/delta_t^2) - (C2/(2*delta_t));
b2     = KK(2) - (2*MM(2)/(delta_t^2));
p2     = data(:)' * MM(2);
p_hat2 = zeros(1, numel(tt)); p_hat2(1) = p2(1);
U2(2)  = (p_hat2(1)/K_hat2)*scalePGA;
for i=2:400
    p_hat2(i) = p2(i) - a2*U2(i-1) - b2*U2(i);
    U2(i+1)   = (p_hat2(i)/K_hat2)*scalePGA;
end

% Mode 3
C3     = 2*Wn(3)*MM(3)*zeta_i(3);
U3     = zeros(1, numel(tt)); U3(1) = 0;
K_hat3 = (MM(3)/(delta_t^2)) + (C3/(2*delta_t));
a3     = (MM(3)/delta_t^2) - (C3/(2*delta_t));
b3     = KK(3) - (2*MM(3)/(delta_t^2));
p3     = data(:)' * MM(3);
p_hat3 = zeros(1, numel(tt)); p_hat3(1) = p3(1);
U3(2)  = (p_hat3(1)/K_hat3)*scalePGA;
for i=2:400
    p_hat3(i) = p3(i) - a3*U3(i-1) - b3*U3(i);
    U3(i+1)   = (p_hat3(i)/K_hat3)*scalePGA;
end

% Mode 4
C4     = 2*Wn(4)*MM(4)*zeta_i(4);
U4     = zeros(1, numel(tt)); U4(1) = 0;
K_hat4 = (MM(4)/(delta_t^2)) + (C4/(2*delta_t));
a4     = (MM(4)/delta_t^2) - (C4/(2*delta_t));
b4     = KK(4) - (2*MM(4)/(delta_t^2));
p4     = data(:)' * MM(4);
p_hat4 = zeros(1, numel(tt)); p_hat4(1) = p4(1);
U4(2)  = (p_hat4(1)/K_hat4)*scalePGA;
for i=2:400
    p_hat4(i) = p4(i) - a4*U4(i-1) - b4*U4(i);
    U4(i+1)   = (p_hat4(i)/K_hat4)*scalePGA;
end

% Mode 5
C5     = 2*Wn(5)*MM(5)*zeta_i(5);
U5     = zeros(1, numel(tt)); U5(1) = 0;
K_hat5 = (MM(5)/(delta_t^2)) + (C5/(2*delta_t));
a5     = (MM(5)/delta_t^2) - (C5/(2*delta_t));
b5     = KK(5) - (2*MM(5)/(delta_t^2));
p5     = data(:)' * MM(5);
p_hat5 = zeros(1, numel(tt)); p_hat5(1) = p5(1);
U5(2)  = (p_hat5(1)/K_hat5)*scalePGA;
for i=2:400
    p_hat5(i) = p5(i) - a5*U5(i-1) - b5*U5(i);
    U5(i+1)   = (p_hat5(i)/K_hat5)*scalePGA;
end

U = [U1;U2;U3;U4;U5];
A = zeros(nDOF, numel(tt)); % modal accelerations
for i=1:nDOF
    A(i,:) = (Wn(i))^2 .* U(i,:);
end

%% ---------------------------
%% Compute responses from modal combination (El-Centro)
Xtt = zeros(nDOF,nDOF);   % shape factors for displacement reconstruction
for i=1:nDOF
    Xtt(:,i) = Gama(i) .* phi(:,i);
end

% Story displacement time histories via modal superposition:
X = zeros(nDOF, numel(tt));
for i=1:nDOF
    % X(i,:) = sum_j Xtt(i,j)*U(j,:)
    X(i,:) = Xtt(i,1).*U(1,:) + Xtt(i,2).*U(2,:) + Xtt(i,3).*U(3,:) ...
           + Xtt(i,4).*U(4,:) + Xtt(i,5).*U(5,:);
end

Xmax = zeros(1,nDOF);
for i=1:nDOF
    Xmax(i) = max(abs(X(i,:)));
end

% Alternative reconstruction using influence form:
U_st = zeros(nDOF,nDOF);
for i=1:nDOF
    U_st(:,i) = (Gama(i)/Wn(i)^2) .* phi(:,i);
end
Uroof_m1 = U_st(:,1)*A(1,:); % all-story contributions from mode 1 etc.
Uroof_m2 = U_st(:,2)*A(2,:);
Uroof_m3 = U_st(:,3)*A(3,:);
Uroof_m4 = U_st(:,4)*A(4,:);
Uroof_m5 = U_st(:,5)*A(5,:);

% Base shear & OTM from modal accelerations
Vbbb = zeros(nDOF, numel(tt));
for i=1:nDOF
    Vbbb(i,:) = Mstr(i) * A(i,:); % modal base shear time history
end
VB    = sum(Vbbb, 1);
VBmax = max(abs(VB));

Mbbb = zeros(nDOF, numel(tt));
for i=1:nDOF
    Mbbb(i,:) = Hstr(i) * Mstr(i) * A(i,:);
end
MB    = sum(Mbbb, 1);
MBmax = max(abs(MB));

% Interstory drift influence matrix (user-provided; fixed typos)
phidelta = [ 0.0146 -0.0392  0.0514  0.0472 -0.0281;
             0.0135 -0.0122 -0.0386  0.0864  0.0753;
             0.0111  0.0233  0.0618  0.0246  0.0986;
             0.0800  0.0427  0.0191  0.0660  0.0906;
             0.0420  0.0326  0.0673 -0.0795  0.0538 ];

drift_st = zeros(nDOF,nDOF);
for i=1:nDOF
    drift_st(i,:) = (Gama(i)/Wn(i)^2) .* phidelta(i,:);
end
Droof_m1 = drift_st(:,1)*A(1,:);
Droof_m2 = drift_st(:,2)*A(2,:);
Droof_m3 = drift_st(:,3)*A(3,:);
Droof_m4 = drift_st(:,4)*A(4,:);
Droof_m5 = drift_st(:,5)*A(5,:);
Droof    = Droof_m1 + Droof_m2 + Droof_m3 + Droof_m4 + Droof_m5;

Droofmax = zeros(1,nDOF);
for i=1:nDOF
    Droofmax(i) = max(abs(Droof(i,:)));
end

disp('****** Answers of Problem 5 (El-Centro time history) ******');
disp('Roof displacement X(5,:):'); disp(X(5,:));
disp('Roof drift Droof(5,:):');   disp(Droof(5,:));
disp('Base shear VB (time history):'); disp(VB);
disp('OTM MB (time history):');       disp(MB);
disp('Base shear envelope (|VB|_max):'); disp(VBmax);
disp('OTM envelope (|MB|_max):');       disp(MBmax);

figure(4)
subplot(2,3,1); plot(tt, X(5,:));    grid on; xlabel('t [s]'); ylabel('u_{roof} [in]'); title('Roof Disp (El-Centro)');
subplot(2,3,2); plot(tt, Droof(5,:));grid on; xlabel('t [s]'); ylabel('\Delta_{roof} [in]'); title('Roof Drift (El-Centro)');
subplot(2,3,3); plot(tt, VB);        grid on; xlabel('t [s]'); ylabel('V_b [kips]'); title('Base Shear (El-Centro)');
subplot(2,3,4); plot(tt, MB);        grid on; xlabel('t [s]'); ylabel('M_b [kip-ft]'); title('OTM (El-Centro)');
subplot(2,3,5); plot(abs(Xmax), H);  grid on; xlabel('|u_{max}| [in]'); ylabel('Story'); title('Max Story Disp (El-Centro)');
subplot(2,3,6); plot(abs(Droofmax),H); grid on; xlabel('|drift_{max}| [in]'); ylabel('Story'); title('Max Story Drift (El-Centro)');

%% ================================================================
%% Response Spectrum–style envelopes using max |U| and |A|
%% ================================================================
Umax = zeros(1,nDOF);
Amax = zeros(1,nDOF);
for i=1:nDOF
    Umax(i) = max(abs(U(i,:)));
    Amax(i) = max(abs(A(i,:)));
end

% Displacement envelope by summing modal contributions of max values
Xx = zeros(1,nDOF);
for i=1:nDOF
    Xx(i) = Xtt(i,1)*Umax(1) + Xtt(i,2)*Umax(2) + Xtt(i,3)*Umax(3) ...
          + Xtt(i,4)*Umax(4) + Xtt(i,5)*Umax(5);
end

% Base shear and OTM envelopes from peak modal accelerations
Vbbbb = zeros(1,nDOF);
Mbbbb = zeros(1,nDOF);
for i=1:nDOF
    Vbbbb(i) = Mstr(i)*Amax(i);
    Mbbbb(i) = Hstr(i)*Mstr(i)*Amax(i);
end
VBB = sum(Vbbbb);
MBB = sum(Mbbbb);

Ddroof_m1 = drift_st(:,1)*Amax(1);
Ddroof_m2 = drift_st(:,2)*Amax(2);
Ddroof_m3 = drift_st(:,3)*Amax(3);
Ddroof_m4 = drift_st(:,4)*Amax(4);
Ddroof_m5 = drift_st(:,5)*Amax(5);
Drooof    = Ddroof_m1 + Ddroof_m2 + Ddroof_m3 + Ddroof_m4 + Ddroof_m5;

disp('****** Answers of Problem 6 (Envelope) ******');
disp('Base shear envelope VBB:'); disp(VBB);
disp('OTM envelope MBB:');       disp(MBB);
disp('Story drift envelope (vector Drooof):'); disp(Drooof);

figure(5)
subplot(1,2,1); plot(abs(Xx), H);     grid on; xlabel('|u_{env}| [in]'); ylabel('Story'); title('Max Story Disp (Envelope)');
subplot(1,2,2); plot(abs(Drooof), H); grid on; xlabel('|drift_{env}| [in]'); ylabel('Story'); title('Max Story Drift (Envelope)');

%% ---------------------------
%% Problem 7 (reporting modal effective mass & height)
disp('****** Answers of Problem 7 ******');
disp('Modal base shear factors (Mstr):'); disp(Mstr);
disp('Effective modal heights (Hstr):');  disp(Hstr);