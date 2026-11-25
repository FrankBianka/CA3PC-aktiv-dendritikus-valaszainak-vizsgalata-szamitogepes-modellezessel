# Szakdolgozat melléklet – Hippokampális CA3 piramissejt MATLAB szimuláció

Ez a repository a hippokampális CA3 piramissejtek dendritikus kalcium-spike-jait vizsgáló szakdolgozat szimulációs kódjait és segédfüggvényeit tartalmazza.

## Modell és módszertan röviden
- **Kiindulási modell:** Kétkompartmentes Pinsky–Rinzel modell (Paul Miller-féle megvalósítás alapján).
- **Saját módosítások és kiegészítések:**
  - `PR_Na0_aramlepcso.m`: az alapmodell kiegészítése $\text{Na}^+$-blokkolással ($G_{\text{Na}} = 0$) és dendritikus áramlépcsővel.
  - `pr_cal_ghk_simulation.m` és `gKCa_shelldepth.m`: dendritikus L-típusú $\text{Ca}^{2+}$-csatorna aktivációs kinetika és Goldman–Hodgkin–Katz (GHK) áramegyenlet, szubmembrán kalciumdinamika (`shell_depth`), valamint $G_{\text{KCa}}$ skálázás és paraméterrácsvizsgálat.
- **Alapértelmezett szimulációs beállítások:** TTX-analóg ($G_{\text{Na}} = 0$), időlépés $dt = 1\,\mu\text{s}$, szimulációs idő $t_{\max} = 2\text{ s}$, dendritikus lépcsőzetes áraminjekció ($0.5\text{–}1.5\text{ s}$).

## Repository felépítése és fájlok

### Fő szimulációs szkriptek
- `src/PR_Na0_aramlepcso.m` – A Paul Miller-féle alapmodellre épülő, $\text{Na}^+$-blokkolással ($G_{\text{Na}} = 0$) és dendritikus áramlépcsővel kiegészített szimulációs szkript a küszöbértékek és válaszok vizsgálatához.
- `src/pr_cal_ghk_simulation.m` – Fő szimulációs szkript a szomatikus és dendritikus membránpotenciálok ($V_S, V_D$), az intracelluláris kalciumkoncentráció ($[\text{Ca}^{2+}]$) és a dendritikus ionáramok ($I_{\text{CaL}}, I_{\text{KCa}}, I_{\text{KAHP}}$) GHK-alapú számításához és ábrázolásához.
- `src/gKCa_shelldepth.m` – Kétdimenziós paraméterrácsvizsgálat ($G_{\text{KCa}} \times \text{shell\_depth}$). TTX-analóg ($G_{\text{Na}} = 0$) körülmények között, $0.8\,\text{nA}$-es dendritikus áraminjekcióval ($0.5\text{–}1.5\,\text{s}$) 30 paraméterkombinációt futtat le (`kca_mults = [0.25, 0.5, 1, 2, 4]`, `shell_mults = [1, 2, 5, 10, 20, 50]`), majd kirajzolja a maximális dendritikus membránpotenciál és a maximális $[\text{Ca}^{2+}]$ hőtérképeit, valamint kiírja a `T_summary` összegző táblázatot.

### Segédfüggvények (kapuzási kinetika)
- `src/PR_dend_gating.m` – A dendritikus kompartment ioncsatornáinak feszültségfüggő kapuzási sebességi állandóit és egyensúlyi értékeit számító segédfüggvény.
- `src/PR_soma_gating.m` – A szomatikus kompartment ioncsatornáinak kapuzási kinetikáját számító segédfüggvény.

## Követelmények és használat
1. **Környezet:** MATLAB (R2021a vagy újabb ajánlott).
2. **Futtatás:** 
  - A repository klónozható, illetve a fájlok külön is letölthetők.
  - A `src/` mappa megnyitandó a MATLAB *Current Folder* ablakában (hogy a segédfüggvények elérhetők legyenek).
  - A kívánt szkript (pl. `PR_Na0_aramlepcso.m` vagy `pr_cal_ghk_simulation.m`) megnyitása után a paraméterek (pl. `amplitudes`, `P_CaL`, `G_KCa`, `shell_depth`) beállíthatók, majd a szimuláció a **Run** gombbal indítható el.
