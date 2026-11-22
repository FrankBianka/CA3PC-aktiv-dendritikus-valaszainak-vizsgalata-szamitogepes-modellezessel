clear

% Időlépés és idővektor definíciója
dt = 1e-6; %10e-6
tmax = 2;
N = round(tmax/dt);
t = (0:N-1)*dt;


% Reversal potenciálok és fizikai konstansok
E_L = -0.060;   % leak reversal potential
E_Na = 0.060;   % sodium
E_K = -0.075;   % potassium
E_Ca = 0.080;   % calcium
F = 96485;  
Rg = 8.314;  
T = 273.15 + 34;      
Ca_out  = 2e-3;  
Ca_rest = 5e-8;

% Kapacitások és felületek
c_m    = 1e-6; % F/cm^2
Cm_tot = 100e-12; % F

S_frac = 1/3;  % soma részarány
D_frac = 1 - S_frac; % dendrit részarány

CmS = Cm_tot * S_frac;
CmD = Cm_tot * D_frac;

area_D = Cm_tot / c_m * D_frac; % cm^2

% Kalcium térfogat és átváltási tényező
shell_depth = 2e-6;                     % cm eredeti 2e-6
vol_L       = area_D * shell_depth * 1e-3; % liter
convert_Ca  = 1 / (2 * F * vol_L);

% Vezetőképességek
G_LS = 5e-9 * S_frac;     
G_Na = 0; %3e-6 * S_frac     
G_K  = 2e-6 * S_frac;      

G_LD   = 5e-9 * D_frac;         
G_Ca   = 2e-6 * D_frac;         
G_KAHP = 0.04e-6 * D_frac;    
G_KCa  = 2.5e-6 * D_frac;      
G_Link = 20e-9; %20e-9

% L-típusú permeabilitás
P_CaL = 0.001;
tau_Ca = 50e-3;             

% Tömbök inicializálása
[VS, VD, Ca, m, h, n, mL, mkca, mkahp] = deal(zeros(1, N));
[I_CaL, I_KCa, I_KAHP, I_Link]          = deal(zeros(1, N));
[I_LS, I_LD, I_Na, I_K, IS, ID] = deal(zeros(1,N));

% Kezdeti feltételek
VS(1) = E_L;  
VD(1) = E_L;  
Ca(1) = Ca_rest; 

[am, bm, ah, bh, an, bn] = soma_gating(E_L);
m(1) = am/(am+bm);  h(1) = ah/(ah+bh);  n(1) = an/(an+bn);

[akca, bkca, akahp, bkahp] = dend_gating(E_L, Ca(1));
mkca(1)  = akca/(akca+bkca);
mkahp(1) = akahp/(akahp+bkahp);

[aL, bL] = CaL_gating(E_L);
mL(1) = aL/(aL+bL);

upd = @(x, a, b) a/(a+b) + (x - a/(a+b))*exp(-dt*(a+b));


%Lépcsőzetes áraminjekció
I_inj_dend = zeros(1, N);
t_start = 0.5;         % Kezdő időpont (sec)
t_end = 1.5;           % Záró időpont (sec)
amplitudes = [0.8e-9]; % Lépcső amplitúdója (nA)
num_steps = length(amplitudes);
step_dur = (t_end - t_start) / num_steps;       % Egy lépcső hossza

for i = 1:N
    curr_t = t(i);
    if curr_t >= t_start && curr_t < t_end
        % Kiszámítjuk, hogy hanyadik lépcsőben vagyunk
        step_idx = floor((curr_t - t_start) / step_dur) + 1;
        if step_idx > num_steps
            step_idx = num_steps;
        end
        I_inj_dend(i) = amplitudes(step_idx);
    end
end

% Szimulációs ciklus
for i = 2:N
    Vs  = VS(i-1);  
    Vd  = VD(i-1);  
    Cai = Ca(i-1);
    
    [am, bm, ah, bh, an, bn] = soma_gating(Vs);
    [akca, bkca, akahp, bkahp] = dend_gating(Vd, Cai);
    [aL, bL]                 = CaL_gating(Vd);
   
    m(i)     = upd(m(i-1), am, bm);
    h(i)     = upd(h(i-1), ah, bh);
    n(i)     = upd(n(i-1), an, bn);
    mL(i)    = upd(mL(i-1), aL, bL);
    mkca(i)  = upd(mkca(i-1), akca, bkca);
    mkahp(i) = upd(mkahp(i-1), akahp, bkahp);
     
    G_Na_now = G_Na * m(i)^2 * h(i);
    G_K_now  = G_K * n(i)^2;
   
    % L-típusú Ca áram (GHK)
    Cai_mol_cm3 = max(Cai, 0) * 1e-3;
    Cao_mol_cm3 = Ca_out * 1e-3;
    J = ghk_Ca(Vd, Cai_mol_cm3, Cao_mol_cm3, T);
    I_CaL(i) = -P_CaL * area_D * (mL(i)^2) * J;
    
    G_KCa_now  = G_KCa * mkca(i) * min(Cai/250e-6, 1);
    G_KAHP_now = G_KAHP * mkahp(i);
    I_KCa(i)   = G_KCa_now * (E_K - Vd);
    I_KAHP(i)  = G_KAHP_now * (E_K - Vd);
    I_Link(i)  = G_Link * (Vd - Vs);
    
    % Leak áramok
    I_LS(i) = G_LS*(E_L - Vs);
    I_LD(i) = G_LD*(E_L - Vd);

    % Szomatikus áramok
    I_Na(i) = G_Na_now*(E_Na - Vs);
    I_K(i)  = G_K_now*(E_K - Vs);

    % Összáramok
    IS(i) = I_LS(i) + I_Na(i) + I_K(i) + I_Link(i);
    ID(i) = I_LD(i) + I_CaL(i) + I_KCa(i) + I_KAHP(i) - I_Link(i) + I_inj_dend(i);

    % Frissítés
    VS(i) = VS(i-1) + dt*IS(i)/CmS;
    VD(i) = VD(i-1) + dt*ID(i)/CmD;
    
    % Kalcium-dinamika frissítés
    Ca_inf = Ca_rest + tau_Ca*convert_Ca*I_CaL(i);
    Ca(i)  = Ca_inf - (Ca_inf - Cai)*exp(-dt/tau_Ca);
end


doPlot = true;
spk = find(VS(1:end-1) < 0 & VS(2:end) >= 0);

if doPlot
    set(0, 'DefaultLineLineWidth', 1.8, ...
           'DefaultLineMarkerSize', 8, ...
           'DefaultAxesLineWidth', 1.5, ...
           'DefaultAxesFontSize', 12, ...
           'DefaultAxesFontWeight', 'Bold');
    
    t_inj_start = 0.45; 
    t_inj_end   = 1.55;
    
    if ~isempty(spk)
        t_spk_start = max(0, t(spk(1)) - 0.01);
        t_spk_end   = min(tmax, t_spk_start + 0.05);
        t_ca_end    = min(tmax, t_spk_start + 0.3);
    else
        t_spk_start = 0.49; t_spk_end = 0.54; t_ca_end = 0.8;
    end

    kca_ymin = min(-0.05, min(I_KCa*1e9)*1.1);
    ca_ymax  = max(1e-6, max(Ca)*1e3*1.1);

    %% FIGURE 1: Teljes időtartam (V_S, V_D)
    figure(1); clf;
    subplot('Position', [0.12 0.58 0.83 0.36])
    plot(t, VS*1000, 'k')
    axis([0 tmax -85 50]); grid on;
    ylabel('V_S (mV)');
    
    subplot('Position', [0.12 0.12 0.83 0.36])
    plot(t, VD*1000, 'k')
    axis([0 tmax -85 50]); grid on;
    ylabel('V_D (mV)'); xlabel('Idő (sec)');

   

  
    %% FIGURE 2: Ionáramok (Kettős Y-tengellyel, közös nullponttal, fekete tengelyekkel)
figure(2); clf;
ax = gca;
xlim([t_inj_start t_inj_end]);
idx_win = (t >= t_inj_start) & (t <= t_inj_end);

% Bal oldali Y-tengely (I_CaL, I_KCa)
yyaxis left
plot(t, I_CaL * 1e9, 'r', 'LineWidth', 1.5); hold on;
plot(t, I_KCa * 1e9, 'g', 'LineWidth', 1.5);
ylabel('I_{CaL}, I_{KCa} (nA)');
ax.YColor = 'k';

yL_vals = [I_CaL(idx_win)*1e9, I_KCa(idx_win)*1e9];
yL_min = min(0, min(yL_vals)*1.1);
yL_max = max(0, max(yL_vals)*1.1);
if yL_min == yL_max, yL_min = -1; yL_max = 1; end
ylim([yL_min, yL_max]);

% Jobb oldali Y-tengely (I_KAHP)
yyaxis right
plot(t, I_KAHP * 1e12, 'b--', 'LineWidth', 1.5);
ylabel('I_{KAHP} (pA)');
ax.YColor = 'k';

yR_vals = I_KAHP(idx_win)*1e12;
yR_min_raw = min(0, min(yR_vals)*1.1);
yR_max_raw = max(0, max(yR_vals)*1.1);
if yR_min_raw == yR_max_raw, yR_min_raw = -1; yR_max_raw = 1; end

% A 0 pozíciójának illesztése a bal tengely arányához:
% zero_frac = (0 - yL_min) / (yL_max - yL_min)
zero_frac = (0 - yL_min) / (yL_max - yL_min);

% Olyan [yR_min, yR_max] tartományt keresünk, amely lefedi a raw határokat
% és amelyben -yR_min / (yR_max - yR_min) == zero_frac
span_from_max = (yR_max_raw > 0) * (yR_max_raw / (1 - zero_frac));
span_from_min = (yR_min_raw < 0) * (-yR_min_raw / zero_frac);
yR_span = max([span_from_max, span_from_min, 1e-6]);

yR_min = -zero_frac * yR_span;
yR_max = (1 - zero_frac) * yR_span;
ylim([yR_min, yR_max]);

grid on; xlabel('Idő (sec)');
title('Ionáramok dinamikája (Kettős Y-tengely)');
legend('I_{CaL} (L-típusú Ca)', 'I_{KCa} (KCa)', 'I_{KAHP} (KAHP)', 'Location', 'best');

    %% FIGURE 3: Külön alablakos Ionáramok
    figure(3); clf;
    subplot(3,1,1)
    plot(t, I_CaL * 1e9, 'r')
    ylabel('I_{CaL} (nA)'); grid on;
    title('Ionáramok komponensei'); xlim([t_spk_start t_spk_end]);
    
    subplot(3,1,2)
    plot(t, I_KCa * 1e9, 'g')
    ylabel('I_{KCa} (nA)'); grid on;
    xlim([t_spk_start t_spk_end]);
    
    subplot(3,1,3)
    plot(t, I_KAHP * 1e12, 'b')
    ylabel('I_{KAHP} (pA)'); xlabel('Idő (sec)'); grid on;
    xlim([t_spk_start t_spk_end]);

    %% FIGURE 4: Membránpotenciál és kalciumszint közelkép
    figure(4); clf;
    subplot('Position',[0.10 0.18 0.38 0.65])
    plot(t,VS*1e3,'k',t,VD*1e3,':k')
    axis([t_spk_start t_spk_end -85 50]);
    grid on;
    ylabel('Membrán Potenciál (mV)');
    xlabel('Idő (sec)');
    legend('V_S','V_D','Location','northeast');

    subplot('Position',[0.58 0.18 0.38 0.65])
    plot(t,Ca*1e3,'k')
    axis([t_spk_start t_spk_end 0 ca_ymax]);
    grid on;
    ylabel('[Ca] (mM)');
    xlabel('Idő (sec)');
    %% FIGURE 5: Membránpotenciál és kalciumszint zoomolt közelkép (0.4 - 0.6 s)
figure(5); clf;
t_zoom_start = 0.4;
t_zoom_end   = 0.6;
idx_zoom = (t >= t_zoom_start) & (t <= t_zoom_end);

% Bal oldali subplot: Membránpotenciál
subplot('Position',[0.10 0.18 0.38 0.65])
plot(t, VS*1e3, 'k', t, VD*1e3, ':k')
xlim([t_zoom_start t_zoom_end]);
ylim([-85 50]);
xticks([0.4 0.5 0.6]);
grid on;
ylabel('Membrán Potenciál (mV)');
xlabel('Idő (sec)');
legend('V_S', 'V_D', 'Location', 'northeast');

% Jobb oldali subplot: Kalciumszint az adott időablakra skálázva
subplot('Position',[0.58 0.18 0.38 0.65])
plot(t, Ca*1e3, 'k')
xlim([t_zoom_start t_zoom_end]);
ca_zoom_max = max(Ca(idx_zoom)*1e3);
if isempty(ca_zoom_max) || ca_zoom_max <= 0
    ca_zoom_max = 1e-3;
else
    ca_zoom_max = ca_zoom_max * 1.15;
end
ylim([0 ca_zoom_max]);
xticks([0.4 0.5 0.6]);
grid on;
ylabel('[Ca] (mM)');
xlabel('Idő (sec)');
end

%%Segédfüggvények
function J = ghk_Ca(V, Ci, Co, T)
zF = 2*96485;
u  = zF*V/(8.314*T);
if abs(u) < 1e-6
    J = zF*(Ci - Co);
else
    J = zF*u*(Ci - Co*exp(-u))/(-expm1(-u));
end
end

function [aL, bL] = CaL_gating(Vd)
v  = Vd*1e3;
aL = 300*exp((v - 6.5)/17);
bL = 300*exp(-(v - 6.5)/17);
end

function [am, bm, ah, bh, an, bn] = soma_gating(Vm)
x = Vm + 0.0469;
if abs(x) < 1e-9, am = 1280; else, am = 320e3*x/(1 - exp(-250*x)); end
x = Vm + 0.0199;
if abs(x) < 1e-9, bm = 1400; else, bm = 280e3*x/(exp(200*x) - 1); end
ah = 128*exp(-(Vm + 0.043)/0.018);
bh = 4e3/(1 + exp(-200*(Vm + 0.020)));
x = Vm + 0.0249;
if abs(x) < 1e-9, an = 80;   else, an = 16e3*x/(1 - exp(-200*x)); end
bn = 250*exp(-25*(Vm + 0.040));
end

function [akca, bkca, akahp, bkahp] = dend_gating(Vd, Ca)
if Vd > -0.010
    akca = 2e3*exp(-(0.0535 + Vd)/0.027);
    bkca = 0;
else
    akca = exp((Vd + 0.050)/0.011 - (Vd + 0.0535)/0.027)/0.018975;
    bkca = 2e3*exp(-(0.0535 + Vd)/0.027) - akca;
end
akahp = min(20, 20e3*Ca);
bkahp = 4;
end

peakVD_mV   = max(VD)*1000;
peakCa_mM   = max(Ca)*1000;
peakICaL_nA = max(I_CaL)*1e9;
peakIKCa_nA = max(abs(I_KCa))*1e9;


fprintf('\n=== Kvantitatív paraméterek ===\n');
fprintf('Peak VD   = %.2f mV\n',peakVD_mV);
fprintf('Peak Ca   = %.4f mM\n',peakCa_mM);
fprintf('Peak ICaL = %.2f nA\n',peakICaL_nA);
fprintf('Peak IKCa = %.2f nA\n',peakIKCa_nA);

fprintf('VD max: %.2f mV, VD vege: %.2f mV\n', max(VD)*1e3, VD(end)*1e3);
fprintf('ICaL min: %.2f nA, ICaL max: %.2f nA\n', min(I_CaL)*1e9, max(I_CaL)*1e9);

peakIinj_nA = max(I_inj_dend) * 1e9;
if any(I_inj_dend > 0)
    inj_amp_nA = mean(I_inj_dend(I_inj_dend > 0)) * 1e9;
else
    inj_amp_nA = 0;
end

fprintf('Peak I_inj    = %.2f nA (amplitudo: %.2f nA, %.2f - %.2f s)\n', ...
    peakIinj_nA, inj_amp_nA, t_start, t_end);