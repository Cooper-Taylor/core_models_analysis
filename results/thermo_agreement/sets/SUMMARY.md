# Three-source ΔG agreement: named reaction sets

Built from `results/thermo_agreement/reaction_features.tsv` (8,359 reactions carrying all three sources' ΔG′°).

## `set_A_three_way_concordant.tsv` — 2,722 reactions

all three sources within 5 kcal/mol and at least one |ΔG′°| > 2 kcal/mol (trivial isomerase zeros excluded), KEGG mapping vouched; sorted largest-ΔG first

| reaction | name | GC | eQ | dGP |
|---|---|---:|---:|---:|
| rxn16706 | 9beta-pimara-7,15-diene,NADPH:oxygen 19-oxidoreductase | -318.40 | -316.94 | -315.67 |
| rxn37184 | methylsterol monooxygenase | -314.72 | -316.94 | -315.09 |
| rxn40150 | 11-oxo-beta-amyrin,NADPH:oxygen oxidoreductase (30-hydroxy | -314.72 | -316.94 | -313.03 |
| rxn11683 | R07509 | -314.72 | -316.94 | -313.33 |
| rxn37185 | methylsterol monooxygenase | -314.72 | -316.93 | -316.84 |
| rxn14236 | 7alpha,26-dihydroxy-4-cholesten-3-one,NADPH:oxygen oxidore | -222.14 | -222.82 | -218.77 |
| rxn14371 | cholest-5-ene-3beta,26-diol,NADPH:oxygen oxidoreductase (2 | -222.14 | -222.82 | -218.77 |
| rxn39787 | Trichloroethanol, NADPH:oxygen oxidoreductase (RH-hydroxyl | -222.14 | -222.81 | -218.71 |

## `set_B_group_transfer_discordant.tsv` — 385 reactions

group-transfer chemistry, vouched mapping, sources spread >20 kcal/mol

| reaction | name | GC | eQ | dGP |
|---|---|---:|---:|---:|
| rxn04474 | R06625 | 1706.37 | 1092.32 | -1.14 |
| rxn07706 | geranylgeranyl-diphosphate:geranylgeranyl-diphosphate gera | 585.49 | 873.99 | -37.19 |
| rxn31355 | PSY (PHYTOENE SYNTHASE); geranylgeranyl-diphosphate gerany | 585.49 | 873.99 | -37.19 |
| rxn04305 | R06443 | 555.27 | 414.14 | 33.25 |
| rxn15678 | R07456 | 398.26 | 293.01 | -75.63 |
| rxn04283 | R06421 | 400.30 | 389.30 | -8.05 |
| rxn40257 | ATP:L-threonyl,bicarbonate adenylyltransferase | -17.74 | 289.39 | 26.39 |
| rxn33029 | ATP:propanoate adenyltransferase | 11.13 | -281.99 | 19.01 |

## `set_B2_group_transfer_discordant_ordinary_scale.tsv` — 352 reactions

as set B but restricted to reactions where every source is inside ±100 kcal/mol -- ordinary metabolic chemistry, not aggregate/polymer reactions

| reaction | name | GC | eQ | dGP |
|---|---|---:|---:|---:|
| rxn20834 | GTP cyclohydrolase IV | -21.13 | 63.93 | -44.08 |
| rxn19830 | GLUTAMIDOTRANS-RXN.d | -25.53 | 58.35 | -3.01 |
| rxn03135 | R04558 | 25.53 | -58.35 | -3.01 |
| rxn02640 | Choloyl-CoA:glycine N-choloyltransferase | -5.29 | -6.64 | -86.81 |
| rxn25846 | RXN-9800.x | -5.29 | -6.64 | -86.81 |
| rxn10475 | R08217 | -14.03 | -2.08 | 66.15 |
| rxn07777 | malonyl-CoA:cinnamoyl-CoA malonyltransferase (cyclizing) | 8.92 | -42.95 | 35.10 |
| rxn07778 | malonyl-CoA:caffeoyl-CoA malonyltransferase (cyclizing) | 8.92 | -42.88 | 35.12 |

## `set_C_kegg_mismapped_withheld.tsv` — 17,271 reactions

the reactions the mask withholds: their stored dGPredictor ΔG′° was predicted from a KEGG reaction ModelSEED does not list for them. Sorted by how many ModelSEED reactions share that KEGG id. `dg_dgp` is the withheld value; `dg_gc`/`dg_eq` are blank where the reaction has no Group-Contribution / eQuilibrator value either.

| reaction | name | GC | eQ | dGP |
|---|---|---:|---:|---:|
| rxn13478 | Heptadecanoate-transport-via-proton-symport | nan | nan | -67.33 |
| rxn13479 | Isomerase-for-keto-meroacid-2 | nan | nan | -67.33 |
| rxn13480 | Isomerase-for-methoxy-meroacid-2 | nan | nan | -67.33 |
| rxn13481 | keto-ylation-1-for-keto-meroacid-1 | nan | nan | -67.33 |
| rxn13482 | keto-ylation-1-for-keto-meroacid-2 | nan | nan | -67.33 |
| rxn13483 | keto-ylation-2-for-keto-meroacid-1 | nan | nan | -67.33 |
| rxn13484 | keto-ylation-2-for-keto-meroacid-2 | nan | nan | -67.33 |
| rxn13485 | Keto-ylation-1-for-keto-meroacid-1 | nan | nan | -67.33 |

## `set_D_all_three_disagree.tsv` — 1,374 reactions

Group Contribution and eQuilibrator themselves differ by >15 kcal/mol

| reaction | name | GC | eQ | dGP |
|---|---|---:|---:|---:|
| rxn04474 | R06625 | 1706.37 | 1092.32 | -1.14 |
| rxn03910 | 2-Phospho-4-(cytidine 5'-diphospho)-2-C-methyl-D-erythrito | 26.77 | 432.65 | 7.30 |
| rxn28646 | 2-C-methyl-D-erythritol 2,4-cyclodiphosphate synthase | 26.77 | 432.65 | 7.30 |
| rxn40257 | ATP:L-threonyl,bicarbonate adenylyltransferase | -17.74 | 289.39 | 26.39 |
| rxn00986 | ATP:propanoate adenylyltransferase | -11.13 | 281.99 | 19.01 |
| rxn35418 | ATP:propanoate adenyltransferase | -11.13 | 281.99 | 19.01 |
| rxn33029 | ATP:propanoate adenyltransferase | 11.13 | -281.99 | 19.01 |
| rxn31212 | acetyl-CoA synthetase (acetate-CoA ligase), putative | -9.59 | 281.50 | 13.94 |
