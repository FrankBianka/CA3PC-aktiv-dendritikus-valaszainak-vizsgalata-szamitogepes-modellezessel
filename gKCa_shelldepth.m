clear; clc;

%% 1. Alap paraméterek és fizikai konstansok
dt   = 1e-6;            % s
tmax = 2.0;             % s
N    = round(tmax/dt);
t    = (0:N-1)*dt;

E_L  = -0.060;          % V (leak reversal potential)
E_Na =  0.060;          % V
E_K  = -0.075;          % V
E_Ca =  0.080;          % V
F    = 96485;           % C/mol
Rg   = 8.314;           % J/(mol*K)
T    = 273.15 + 34;     % K
Ca_out  = 2e-3;         % M (2 mM)
Ca_rest = 5e-8;         % M (50 nM)

c_m    = 1e-6;          % F/cm^2
Cm_tot = 100e-12;       % F

S_frac = 1/3;           % soma arány
D_frac = 1 - S_frac;    % dendrit arány

CmS    = Cm_tot * S_frac;
CmD    = Cm_tot * D_frac;
area_D = (Cm_tot / c_m) * D_frac; % cm^2

% Alap vezetőképességek
G_LS       = 5e-9 * S_frac;
G_Na       = 0;             % Na blokkolva
G_K        = 2e-6 * S_frac;

G_LD       = 5e-9 * D_frac;
G_Ca       = 2e-6 * D_frac;
G_KAHP     = 0.04e-6 * D_frac;
G_KCa_base = 2.5e-6 * D_frac;
G_Link     = 20e-9;

% L-típusú Ca permeabilitás és Ca-időállandó
P_CaL  = 0.001;
tau_Ca = 50e-3;

% Alap shell_depth (cm)
shell_depth_base = 2e-6; % cm

% Lépcsőzetes dendritikus áraminjekció
I_inj_dend = zeros(1, N);
t_start    = 0.5;       % s
t_end      = 1.5;       % s
inj_amp    = 0.8e-9;    % A
I_inj_dend(t >= t_start & t < t_end) = inj_amp;

%% 2. Paraméterrács definíciója
kca_mults   = [0.25, 0.5, 1, 2, 4];
shell_mults = [1, 2, 5, 10, 20, 50];

nKCa   = numel(kca_mults);
nShell = numel(shell_mults);

% Eredménymátrixok: sorok = G_KCa szorzók, oszlopok = shell_depth szorzók
max_VD_mV   = nan(nKCa, nShell);   % abszolút max V_D (mV)
delta_VD_mV = nan(nKCa, nShell);   % max(V_D, injekció alatt) - V_D(t_start-dt) (mV)
max_Ca_mM   = nan(nKCa, nShell);   % mM
max_Ca_uM   = nan(nKCa, nShell);   % uM

% Nyugalmi gating értékek kiszámítása
[am0, bm0, ah0, bh0, an0, bn0] = soma_gating(E_L);
m0 = am0/(am0+bm0); h0 = ah0/(ah0+bh0); n0 = an0/(an0+bn0);
[aL0, bL0] = CaL_gating(E_L);
mL0 = aL0/(aL0+bL0);

upd    = @(x, a, b) a/(a+b) + (x - a/(a+b))*exp(-dt*(a+b));
exp_Ca = exp(-dt/tau_Ca);

idx_inj = (t >= t_start) & (t < t_end);
idx_pre = find(t < t_start, 1, 'last');

%% 3. Rács-szimuláció (minden kombinációra)
fprintf('Parameterracs szimulacio inditasa (%d kombinacio)...\n', nKCa * nShell);

for ik = 1:nKCa
    G_KCa_curr = G_KCa_base * kca_mults(ik);
    
    for is = 1:nShell
        shell_depth = shell_depth_base * shell_mults(is);
        vol_L       = area_D * shell_depth * 1e-3; % L
        convert_Ca  = 1 / (2 * F * vol_L);
        
        % Állapotváltozók inicializálása
        VS    = zeros(1, N);
        VD    = zeros(1, N);
        Ca    = zeros(1, N);
        I_CaL = zeros(1, N);
        
        VS(1) = E_L;
        VD(1) = E_L;
        Ca(1) = Ca_rest;
        
        m = m0; h = h0; n = n0; mL = mL0;
        [akca, bkca, akahp, bkahp] = dend_gating(E_L, Ca(1));
        mkca  = akca/(akca+bkca);
        mkahp = akahp/(akahp+bkahp);
        
        % Időléptetéses szimuláció
        for i = 2:N
            Vs  = VS(i-1);
            Vd  = VD(i-1);
            Cai = Ca(i-1);
            
            [am, bm, ah, bh, an, bn]     = soma_gating(Vs);
            [akca, bkca, akahp, bkahp] = dend_gating(Vd, Cai);
            [aL, bL]                     = CaL_gating(Vd);
            
            m     = upd(m, am, bm);
            h     = upd(h, ah, bh);
            n     = upd(n, an, bn);
            mL    = upd(mL, aL, bL);
            mkca  = upd(mkca, akca, bkca);
            mkahp = upd(mkahp, akahp, bkahp);
            
            G_Na_now = G_Na * m^2 * h;
            G_K_now  = G_K * n^2;
            
            Cai_mol_cm3 = max(Cai, 0) * 1e-3;
            Cao_mol_cm3 = Ca_out * 1e-3;
            J = ghk_Ca(Vd, Cai_mol_cm3, Cao_mol_cm3, T);
            I_CaL_now = -P_CaL * area_D * (mL^2) * J;
            I_CaL(i)  = I_CaL_now;
            
            G_KCa_now  = G_KCa_curr * mkca * min(Cai/250e-6, 1);
            G_KAHP_now = G_KAHP * mkahp;
            
            I_KCa_now  = G_KCa_now  * (E_K - Vd);
            I_KAHP_now = G_KAHP_now * (E_K - Vd);
            I_Link_now = G_Link * (Vd - Vs);
            
            I_LS_now = G_LS * (E_L - Vs);
            I_LD_now = G_LD * (E_L - Vd);
            I_Na_now = G_Na_now * (E_Na - Vs);
            I_K_now  = G_K_now  * (E_K - Vs);
            
            IS_now = I_LS_now + I_Na_now + I_K_now + I_Link_now;
            ID_now = I_LD_now + I_CaL_now + I_KCa_now + I_KAHP_now - I_Link_now + I_inj_dend(i);
            
            VS(i) = Vs + dt * IS_now / CmS;
            VD(i) = Vd + dt * ID_now / CmD;
            
            Ca_inf = Ca_rest + tau_Ca * convert_Ca * I_CaL_now;
            Ca(i)  = Ca_inf - (Ca_inf - Cai) * exp_Ca;
        end
        
        %% 4. Metrikák kiértékelése
        vd_pre = VD(idx_pre) * 1e3;
        vd_inj_max = max(VD(idx_inj)) * 1e3;
        
        max_VD_mV(ik, is)   = max(VD) * 1e3;
        delta_VD_mV(ik, is) = vd_inj_max - vd_pre;
        max_Ca_mM(ik, is)   = max(Ca) * 1e3;
        max_Ca_uM(ik, is)   = max(Ca) * 1e6;
        
        fprintf('KCa=%4.2fx | Shell=%2dx -> maxVD=%6.2f mV | dVD=%6.2f mV | maxCa=%.4f mM (%.2f uM)\n', ...
            kca_mults(ik), shell_mults(is), ...
            max_VD_mV(ik, is), delta_VD_mV(ik, is), max_Ca_mM(ik, is), max_Ca_uM(ik, is));
    end
end

%% 5. Hőtérképek: x = shell_depth, y = G_KCa
shell_labels = string(shell_mults) + "x";
kca_labels   = string(kca_mults) + "x";

figure('Name', 'G_KCa es shell_depth hatasa (max V_D es max [Ca])', 'Position', [100 100 1200 480]);

% Bal oldali heatmap: max V_D (mV)
subplot(1, 2, 1);
imagesc(1:nShell, 1:nKCa, max_VD_mV);
set(gca, 'YDir', 'normal', ...
         'XTick', 1:nShell, 'XTickLabel', shell_labels, ...
         'YTick', 1:nKCa,   'YTickLabel', kca_labels, ...
         'FontSize', 11, 'FontWeight', 'bold');
xlabel('Shell depth szorzó (alap: 2 \mum)');
ylabel('G_{KCa} szorzó');
title('Max V_D (mV)');
cb1 = colorbar;
ylabel(cb1, 'Max V_D (mV)');

for ik = 1:nKCa
    for is = 1:nShell
        txt_col = 'w';
        if max_VD_mV(ik, is) < 0.5 * (min(max_VD_mV(:)) + max(max_VD_mV(:)))
            txt_col = 'k';
        end
        text(is, ik, sprintf('%.1f', max_VD_mV(ik, is)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'Color', txt_col, 'FontSize', 9.5, 'FontWeight', 'bold');
    end
end

% Jobb oldali heatmap: max [Ca] (uM)
subplot(1, 2, 2);
imagesc(1:nShell, 1:nKCa, max_Ca_uM);
set(gca, 'YDir', 'normal', ...
         'XTick', 1:nShell, 'XTickLabel', shell_labels, ...
         'YTick', 1:nKCa,   'YTickLabel', kca_labels, ...
         'FontSize', 11, 'FontWeight', 'bold');
xlabel('Shell depth szorzó (alap: 2 \mum)');
ylabel('G_{KCa} szorzó');
title('Max [Ca] (\muM)');
cb2 = colorbar;
ylabel(cb2, 'Max [Ca] (\muM)');

for ik = 1:nKCa
    for is = 1:nShell
        txt_col = 'w';
        if max_Ca_uM(ik, is) < 0.5 * (min(max_Ca_uM(:)) + max(max_Ca_uM(:)))
            txt_col = 'k';
        end
        text(is, ik, sprintf('%.2f', max_Ca_uM(ik, is)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'Color', txt_col, 'FontSize', 9, 'FontWeight', 'bold');
    end
end

%% 6. Összefoglaló táblázat a konzolon
T_summary = table();
row_idx = 1;
for ik = 1:nKCa
    for is = 1:nShell
        T_summary.KCa_mult(row_idx)   = kca_mults(ik);
        T_summary.Shell_mult(row_idx) = shell_mults(is);
        T_summary.Max_VD_mV(row_idx)  = max_VD_mV(ik, is);
        T_summary.dVD_mV(row_idx)     = delta_VD_mV(ik, is);
        T_summary.Max_Ca_uM(row_idx)  = max_Ca_uM(ik, is);
        T_summary.Max_Ca_mM(row_idx)  = max_Ca_mM(ik, is);
        row_idx = row_idx + 1;
    end
end
disp(T_summary);

%% Segédfüggvények
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