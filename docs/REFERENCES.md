# References (PDFs are not committed; obtain them from the sources below)

| Key | Reference | Used for | Copy used |
|---|---|---|---|
| Ayoub2022 | B. Ayoub, S. Moreau, S. Lhostis, H. Frémont, S. Mermoz, E. Souchier, E. Deloffre, S. Escoubas, T. W. Cornelius, O. Thomas, "In-situ characterization of thermomechanical behavior of copper nano-interconnect for 3D integration," *Microelectron. Eng.* 261, 111809 (2022), doi:10.1016/j.mee.2022.111809 | Stage 1 reference: geometry, Cu constants, Ludwick fits, measured ε′zz (Figs. 2–3), Fig. 5 | authors' accepted version, HAL **hal-03672631v1** (https://hal.science/hal-03672631); page numbers in the spec refer to this 7-page file |
| ChangHimmel1966 | Y. A. Chang, L. Himmel, "Temperature dependence of the elastic constants of Cu, Ag, and Au above room temperature," report **UCRL-16697**, Lawrence Radiation Laboratory (1966); published *J. Appl. Phys.* 37, 3567 (1966), doi:10.1063/1.1708903 | C11, C12, C44(T), 300–800 K → G(T) scaling of plasticity (variant 2) | report version (eScholarship); **Table I, report p. 14 (PDF p. 18)** → `inputs/cu_elastic_constants_vs_T.csv` |
| Hahn1970 | T. A. Hahn, "Thermal expansion of copper from 20 to 800 K — Standard Reference Material 736," *J. Appl. Phys.* 41, 5096 (1970), doi:10.1063/1.1658614 | Cu α(T) variant (R09) | **not yet obtained**; `inputs/cu_cte_vs_T.csv.template` |
| LedbetterNaimon1974 | H. M. Ledbetter, E. R. Naimon, "Elastic properties of metals and alloys. II. Copper," *J. Phys. Chem. Ref. Data* 3, 897 (1974), doi:10.1063/1.3253150 | not used: C_ij(T) only as figures | — |
| Nye1985 | J. F. Nye, *Physical Properties of Crystals* (Oxford, 1985) | source of Ayoub's C11 = 168.4, C12 = 121.4, C44 = 75.4 GPa | cited via Ayoub2022 |
