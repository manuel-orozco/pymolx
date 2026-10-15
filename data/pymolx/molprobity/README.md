# MolProbity Top8000 reference data

`top8000.npz` holds the Ramachandran (rama8000) and side-chain rotamer
(rota8000) percentile grids that MolProbity uses for validation, converted
to compressed NumPy arrays by `build_top8000.py`. pymolx uses them in
`structure_diagnostics` (`modules/pymolx/diagnostics.py`): Ramachandran
and rotamer outliers, with MolProbity's cutoffs.

Source: github.com/rlabduke/reference_data, folder `Top8000`, commit
5ee6875fc29eccc3c9dc7cbe705ff9fd1b505d7d (2026-09-14).

Data produced by the Richardson Lab at Duke University, licensed under
the Creative Commons Attribution 4.0 International License (CC BY 4.0,
`LICENSE` here). Copyright (c) 2016-2026 David C. Richardson, Jane S.
Richardson, and contributors. Changes: only the file format (text grids to
float16 arrays).

Please cite:

- Rotamer contour grids: Hintze BJ, Lewis SM, Richardson JS, Richardson
  DC (2016) MolProbity's ultimate rotamer-library distributions for model
  validation. Proteins 84(9):1177-1189. https://doi.org/10.1002/prot.25039
- rama8000 Ramachandran contour grids: Williams CJ, Headd JJ, Moriarty
  NW, Prisant MG, Videau LL, Deis LN, Verma V, Keedy DA, Hintze BJ, Chen
  VB, Jain S, Lewis SM, Arendall WB, Snoeyink J, Adams PD, Lovell SC,
  Richardson JS, Richardson DC (2018) MolProbity: More and better reference
  data for improved all-atom structure validation. Protein Science
  27(1):293-315. https://doi.org/10.1002/pro.3330
