# CA3PC-aktiv-dendritikus-valaszainak-vizsgalata-szamitogepes-modellezessel
A szakdolgozathoz készült MATLAB szimulációs kódok (módosított Pinsky-Rinzel modell, L-típusú Ca2+ áram és GHK-formalizmus)
# Szakdolgozat melléklet – CA3 piramissejt MATLAB szimuláció

Ez a repository a hippokampális CA3 piramissejtek dendritikus kalcium-spike-jait vizsgáló szakdolgozat szimulációs kódját tartalmazza.

## Modell és módszertan röviden
- **Kiindulási modell:** Kétkompartmentes Pinsky–Rinzel modell (Miller-féle variáns).
- **Módosítások:** Dendritikus L-típusú $\text{Ca}^{2+}$-csatorna aktivációs kinetika és Goldman–Hodgkin–Katz (GHK) áramegyenlet, szubmembrán kalciumdinamika (`shell_depth`), valamint $G_{KCa}$ skálázás.
- **Szimulációs beállítás:** TTX-analóg ($G_{Na} = 0$), időlépés $dt = 1\,\mu\text{s}$, $t_{\max} = 2\text{ s}$, dendritikus lépcsőzetes áraminjekció ($0.5\text{–}1.5\text{ s}$).

## Fájlok
- `pr_cal_ghk_simulation.m` – Fő szimulációs szkript a membránpotenciálok ($V_S, V_D$), a kalciumkoncentráció ($[\text{Ca}^{2+}]$) és a dendritikus ionáramok ($I_{CaL}, I_{KCa}, I_{KAHP}$) számításához és ábrázolásához.

## Használat
Fájl megnyitása MATLAB-ban, kívánt paraméterek beállítása (pl. `amplitudes`, `P_CaL`, `G_KCa`, `shell_depth`), majd a szkript futtatása.
