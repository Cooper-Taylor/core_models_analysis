# How each vocabulary encodes the same chemistry

Counts are entries in that vocabulary matching the theme; the examples are verbatim group names.

## phosphate / phosphoanhydride

**fine-tuned dGPredictor** (96 learned fingerprints): `P`, `OP`, `O=P`, `POP`, `COP`, `POS`, `CCP`, `CNP`, `[O-]P`, `O=COP`, `CCCOP`, `O=CCP`, `C=CCOP`, `NC(=O)OP` …

**ModelSEED Group Contribution** (12 groups): `Itriphos`, `PO3`, `WCOOPO3`, `WPO3`, `WPO4nW`, `formylphosphate`, `mid_phos`, `orthophosphate`, `prim_phos`, `pyrophos`, `thiomid_phos`, `thioprim_phos`

**eQuilibrator / component-contribution** (36 groups): `-CO-OPO3 [H0 Z-2 Mg0]`, `-CO-OPO3 [H1 Z-1 Mg0]`, `-CO-OPO3 [H2 Z0 Mg0]`, `CO-OPO3 [H1 Z-2 Mg0]`, `CO-OPO3 [H2 Z-1 Mg0]`, `CO-OPO3 [H3 Z0 Mg0]`, `-OPO3- [H0 Z-1 Mg0]`, `-OPO3- [H1 Z0 Mg0]`, `-OPO2- [H0 Z-1 Mg0]`, `-OPO2- [H1 Z0 Mg0]`, `-OPO3 [H0 Z-2 Mg0]`, `-OPO3 [H1 Z-1 Mg0]`, `-OPO3 [H2 Z0 Mg0]`, `-OPO3-OPO2- [H0 Z-2 Mg0]` …

## carboxylate

**fine-tuned dGPredictor** (182 learned fingerprints): `O=C([O-])O`, `cC(=O)[O-]`, `NC(=O)[O-]`, `CC(=O)[O-]`, `CC(=O)[O-]`, `O=C([O-])O`, `cC(=O)[O-]`, `NC(=O)[O-]`, `O=C([O-])CO`, `CCC(=O)[O-]`, `NCC(=O)[O-]`, `cCC(=O)[O-]`, `O=C([O-])CCO`, `CNCC(=O)[O-]` …

**ModelSEED Group Contribution** (6 groups): `WCOOPO3`, `WCOOn`, `acetate`, `carbamate`, `formate`, `oxalate`

**eQuilibrator / component-contribution** (9 groups): `N-COO [H1 Z0 Mg0]`, `N-COO [H0 Z-1 Mg0]`, `N-COO- [H0 Z0 Mg0]`, `-O-C=O [H0 Z0 Mg0]`, `ring >C=O [H0 Z0 Mg0]`, `-COO [H0 Z-1 Mg0]`, `-COO [H1 Z0 Mg0]`, `>C=O [H0 Z0 Mg0]`, `-C=O [H1 Z0 Mg0]`

## thiol / thioester / disulfide

**fine-tuned dGPredictor** (103 learned fingerprints): `S`, `OS`, `CS`, `O=S`, `POS`, `CSC`, `CCS`, `CSC`, `CCS`, `CSS`, `O=CS`, `O=CS`, `NCCS`, `CCCS` …

**ModelSEED Group Contribution** (26 groups): `FeS`, `S`, `S2-`, `S2O3H`, `S2O4`, `SO2`, `SPb+`, `Sb5+`, `Se`, `Se2-`, `Sn2+`, `Sr2+`, `WSH`, `WSNeg` …

**eQuilibrator / component-contribution** (17 groups): `-S-O [H1 Z0 Mg0]`, `-S-O [H0 Z-1 Mg0]`, `ring -S- [H0 Z0 Mg0]`, `-SO3 [H1 Z0 Mg0]`, `-SO3 [H0 Z-1 Mg0]`, `-SO2- [H0 Z0 Mg0]`, `-SOO [H0 Z-1 Mg0]`, `-SOO [H1 Z0 Mg0]`, `-S< [H0 Z0 Mg0]`, `-S-S- [H0 Z0 Mg0]`, `ring -S-S- [H0 Z0 Mg0]`, `-C(=O)S- [H0 Z0 Mg0]`, `-C(=O)S [H1 Z0 Mg0]`, `-C(=O)S [H0 Z-1 Mg0]` …

## amine / amide

**fine-tuned dGPredictor** (509 learned fingerprints): `N`, `CN`, `NO`, `cN`, `C#N`, `CNC`, `CNP`, `[N]`, `N=O`, `N#N`, `CCN`, `C=N`, `cNC`, `NCN` …

**ModelSEED Group Contribution** (14 groups): `NH4plus`, `RWNHW`, `RWdblNHpW`, `WNH2`, `WNH2W`, `WNH3`, `WNHW`, `WNHWW`, `WdblNH`, `WdblNH2`, `amide`, `hydroxylamine`, `methylamine`, `urea`

**eQuilibrator / component-contribution** (36 groups): `-NO2 [H0 Z0 Mg0]`, `NC(=N)N [H0 Z0 Mg0]`, `two fused rings -N< [H1 Z1 Mg0]`, `two fused rings -N< [H0 Z0 Mg0]`, `two fused rings =N< [H0 Z1 Mg0]`, `ring =N< [H0 Z1 Mg0]`, `-N- [H2 Z1 Mg0]`, `-N- [H1 Z0 Mg0]`, `-N- [H0 Z-1 Mg0]`, `ring =N- [H1 Z1 Mg0]`, `ring =N- [H0 Z0 Mg0]`, `ring -N- [H2 Z1 Mg0]`, `ring -N- [H1 Z0 Mg0]`, `ring -N- [H0 Z-1 Mg0]` …

## aromatic / heteroaromatic

**fine-tuned dGPredictor** (432 learned fingerprints): `cO`, `cF`, `cN`, `cC`, `c=O`, `ncn`, `cnc`, `ccc`, `cCO`, `cCC`, `ccn`, `cCO`, `CCn`, `ncs` …

**ModelSEED Group Contribution** (5 groups): `HeteroAromatic`, `Nap`, `WtrpCH`, `WtrpCW`, `WtrpN`

**eQuilibrator / component-contribution** (60 groups): `ring -Cl [H0 Z0 Mg0]`, `ring -Br [H0 Z0 Mg0]`, `ring -I [H0 Z0 Mg0]`, `ring -F [H0 Z0 Mg0]`, `ring -s- [H0 Z0 Mg0]`, `ring -S- [H0 Z0 Mg0]`, `ring -S-S- [H0 Z0 Mg0]`, `ring -OPO3- [H0 Z-1 Mg0]`, `ring -OPO3- [H1 Z0 Mg0]`, `ring -OPO3-OPO2- [H0 Z-2 Mg0]`, `ring -OPO3-OPO2- [H1 Z-1 Mg0]`, `ring -OPO3-OPO2- [H2 Z0 Mg0]`, `ring -OPO2-OPO2- [H0 Z-2 Mg0]`, `ring -OPO2-OPO2- [H1 Z-1 Mg0]` …

## magnesium

**fine-tuned dGPredictor** (2 learned fingerprints): `[Mg+2]`, `[Mg+2]`

**ModelSEED Group Contribution** (1 groups): `Mg`

**eQuilibrator / component-contribution** (3 groups): `-OPO3-OPO2- [H0 Z0 Mg1]`, `-OPO2-OPO2- [H0 Z0 Mg1]`, `ring -OPO2-OPO2- [H0 Z0 Mg1]`

## Summary counts

| theme | fine-tuned dGPredictor | ModelSEED GC | eQuilibrator |
|---|---:|---:|---:|
| phosphate / phosphoanhydride | 96 | 12 | 36 |
| carboxylate | 182 | 6 | 9 |
| thiol / thioester / disulfide | 103 | 26 | 17 |
| amine / amide | 509 | 14 | 36 |
| aromatic / heteroaromatic | 432 | 5 | 60 |
| magnesium | 2 | 1 | 3 |

## Structural properties of the alphabets

| property | fine-tuned dGPredictor | ModelSEED GC | eQuilibrator |
|---|---|---|---|
| unit | RDKit atom environment, canonical fragment SMILES | named chemical group | named chemical group |
| built by | RDKit, automatically, per compound | MFAToolkit rules | curated SMARTS-style rules |
| protonation state in the label | no — implicit in the pH-7 SMILES (456 learned fragments carry an explicit charge) | no | yes — every group is `[Hn Zq Mgm]` (55 charged variants) |
| Mg binding in the label | no | one free-ion group (`Mg`) | yes (3 Mg-bound variants) |
| ring context in the label | implicit (aromatic lowercase atoms in 432 fragments) | yes (`RW…`, `T…`, `HeteroAromatic`) | yes (54 explicit `ring`/`fused rings` groups) |
| whole-molecule entries | no | yes (`H2O`, `CO2`, `urea`, `acetate`, …) | yes (50 one-hot placeholder columns for non-decomposable compounds) |
| origin / per-molecule constant | none | `Origin` | `Origin [H0 Z0 Mg0]` |
| undecomposable marker | reaction is simply not predicted | `NoGroup` | one-hot placeholder column + RMSE_inf |
