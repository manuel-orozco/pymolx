# Incentive PyMOL parity: gap inventory

This file tracks the work needed for pymolx to look and act as close as
possible to Incentive (commercial) PyMOL. Each gap has an ID, so commits
and PRs can refer to it (for example `parity: P-04 desaturate`).

**Status legend:** `[ ]` todo · `[~]` in progress · `[x]` done · `[-]` won't do (give the reason)

**Effort:** S = up to a day · M = a few days · L = a week or more / touches the C++ core

## Ground rules

- **Clean-room only.** Write the code from public documentation, the
  docstrings already in this repo, and observed behaviour. Don't
  decompile Incentive binaries or copy its icons, splash images or
  bundled Schrödinger libraries.
- **Use our own branding.** "PyMOL" is a Schrödinger trademark. Keep the
  pymolx name, logo and splash (see section H).
- **Keep upstream merges easy.** Put new code in a separate package
  (proposed: `modules/pymolx/`). Replace stubs by rebinding the `cmd`
  names from that package, not by editing upstream files, whenever that
  is possible.
- **Never call an item done without a test.** Many items already have an
  upstream test that is skipped behind `@testing.requires('incentive')`
  (see section I).

## How this list was built

Sections A–D come from markers in the source as of upstream `2a414d80`
(3.2.0a0). They are verified. Sections F–G come from comparing with
Incentive and must be checked side by side before work starts. To
refresh sections A–D after a rebase:

```sh
grep -rn "IncentiveOnlyException\|incentive_format_not_available" modules
grep -rn "_PYMOL_IP_EXTRAS\|_PYMOL_IP_SPLASH\|INCENTIVE_ONLY" layer*
grep -rln "requires('incentive'\|requires(\"incentive\"" testing/tests
```

---

## A. Commands stubbed with `IncentiveOnlyException`

These commands are registered in [keywords.py](modules/pymol/keywords.py)
and documented, but they raise as soon as they are called. The existing
docstrings define the API we have to match.

| ID | Status | Command | Location | Upstream test | Proposed approach | Effort |
|---|---|---|---|---|---|---|
| P-01 | [ ] | `clean` (geometry cleanup / minimisation) | [computing.py:15](modules/pymol/computing.py#L15) | `api/computing.py::testClean`, `jira/PYMOL-2985.py` | RDKit MMFF94 or OpenMM on the selection (honour `present`, `fix`, `restrain`), then write the coordinates back with `load_coords` / `alter_state` | M |
| P-02 | [ ] | `assign_stereo` | [stereochemistry/\_\_init\_\_.py:29](modules/pymol/stereochemistry/__init__.py#L29) | `api/stereochemistry.py` | RDKit CIP labels (`rdCIPLabeler`) stored in the `stereo` atom property. Docstring already names `method=rdkit` | S |
| P-03 | [ ] | `pi_interactions` (pi-pi, pi-cation) | [querying.py:545](modules/pymol/querying.py#L545) | `api/querying.py::testPiInteractions` (1rx1: 4 pi-pi, 1 pi-cation) | Ring perception (RDKit or SSSR on the PyMOL bond graph), centroid and normal geometry, distance and angle cutoffs. Output as a distance object like `distance mode=` | M |
| P-04 | [ ] | `desaturate` | [experimenting.py:280](modules/pymol/experimenting.py#L280) | none | Pure Python: read the colours with `iterate`, blend toward luminance by `a`, write per-atom colours back with `set_color` / `alter` | S |
| P-05 | [ ] | `callout` | [experimenting.py:246](modules/pymol/experimenting.py#L246) | none | Label object plus CGO leader line. Screen-space placement with `screen='auto'`. Check the docstring for the full signature | M |
| P-06 | [ ] | `focal_blur` | [experimenting.py:244](modules/pymol/experimenting.py#L244) | `api/experimenting.py::testFocalblur` | Re-implement in Python: jitter the camera around the focal plane, render N images with `png`, average them with numpy (the classic PyMOL-wiki algorithm by the authors credited in the docstring) | S |
| P-07 | [ ] | `morph` (**all methods**, `linear` too) | [morphing.py:53](modules/pymol/morphing.py#L53) | `api/morphing.py`, `jira/PYMOL-174.py` | Step 1: `method=linear` (sequence match, `align`, then interpolate into states). Step 2: a `rigimol`-like method: rigid-domain detection (DynDom-style), slerp of each domain, then a short `sculpt` refinement (`refinement=`) | S (linear) / L (rigid) |
| P-08 | [ ] | `help_setting` | [helping.py:99](modules/pymol/helping.py#L99) | `api/helping.py::testHelpSetting` | Ship a `settings_docs` data file (name → description, type, default, range) and print from it. Write the text ourselves; the PyMOL Wiki is only a reference (check its licence before copying any text) | M (mostly writing) |
| P-09 | [ ] | `load_mtz` | [importing.py:1511](modules/pymol/importing.py#L1511) | none | `gemmi`: read the MTZ, guess the amplitude/phase columns, FFT to a map, load with `load_brick` / a ccp4 string. Respect `reso_low` / `reso_high` | M |

## B. File formats

| ID | Status | Format | Location | Upstream test | Proposed approach | Effort |
|---|---|---|---|---|---|---|
| F-01 | [ ] | `.mae` / `.maegz` load (incl. `mimic=1` styling, object/atom properties, ANISOU, sub-blocks, group naming) | [importing.py:1620](modules/pymol/importing.py#L1620), [Executive.cpp:3963](layer3/Executive.cpp#L3963) | `api/load_mae.py`, `api/symop.py::test_load_mae`, `api/importing.py::testLoad_mimic`, `properties/*`, `jira/PYMOL-{771,1498,1514,2838,3221}.py` | Python reader built on Schrödinger's open-source **maeparser** (MIT) or a small pure-Python parser, producing a `chempy` model. `mimic=1` maps Maestro display flags to reps and settings | L |
| F-02 | [ ] | `.mae` save | `api/exporting.py:211,271` | same | Writer for `chempy` → MAE. Must handle `chain=" "` (`jira/PYMOL-3231.py`) | M |
| F-03 | [ ] | STL export (`get_stlstr`, `save *.stl`) | [lazyio.py:230](modules/pymol/lazyio.py#L230) | `api/exporting_geom.py::testSTL` | Reuse the triangle data behind the existing COLLADA/OBJ/glTF export paths and write binary STL | S |
| F-04 | [ ] | STL import (`read_stlstr`) | [lazyio.py:240](modules/pymol/lazyio.py#L240) | none | `numpy-stl` or a hand-written parser, producing CGO `TRIANGLES` | S |
| F-05 | [ ] | COLLADA import (`read_collada`) | [lazyio.py:250](modules/pymol/lazyio.py#L250) | none | `pycollada` or `trimesh`, producing CGO | S |
| F-06 | [ ] | `.moe` (MOE molecules/surfaces) | [importing.py:1642](modules/pymol/importing.py#L1642) | `jira/PYMOL-1191.py` | Format is not publicly documented, so it's low priority. Use openly available specs only | L |
| F-07 | [-] | `.vis` (Schrödinger surfaces, HDF5) | [importing.py:1641](modules/pymol/importing.py#L1641) | `jira/PYMOL-1441.py` | Only matters for Maestro users. Revisit after F-01 | L |
| F-08 | [-] | `.phypo` (Phase pharmacophore hypotheses) | [importing.py:1643](modules/pymol/importing.py#L1643) | `api/epymol_ph4.py` | Tied to the Schrödinger Phase product, so skip unless needed | L |
| F-09 | [ ] | `pdb_header` object property on load (`object_props=`) | n/a | `api/importing.py::testPdbHeader` | Confirmed missing: with F-09 unlocked, 6 of 9 `testPdbHeader` cases fail (no `pdb_header` property is attached on load) | S |

## C. C++ core features behind `_PYMOL_IP_EXTRAS`

These are compiled out of the open-source build, and some of the
`#ifdef` blocks are empty. Each one needs a real implementation, not
just a flag change.

| ID | Status | Feature | Location | Upstream test | Notes | Effort |
|---|---|---|---|---|---|---|
| C-01 | [ ] | `sphere_mode` 10 and 11 (higher-quality shader impostor spheres) | [Setting.cpp:2469](layer1/Setting.cpp#L2469) | `cgos/reps.py::test_sphere_mode_10_11` | New sphere shader modes in `data/shaders` + `RepSphere` | L |
| C-02 | [ ] | Volume rendering: `volume_mode` (pre-integrated / alternative modes) | [ObjectVolume.cpp:744](layer2/ObjectVolume.cpp#L744), [:831](layer2/ObjectVolume.cpp#L831), [:945](layer2/ObjectVolume.cpp#L945), [ShaderMgr.cpp:464](layer0/ShaderMgr.cpp#L464) | none | The `volume_mode` setting exists (`SettingInfo.h:849`) but does nothing. Implement pre-integrated transfer-function rendering | L |
| C-03 | [ ] | `stereo` / `text_type` atom properties computed on demand (`iterate ... stereo`, `label`) | [P.cpp:758](layer1/P.cpp#L758), [P.cpp:936](layer1/P.cpp#L936), [ObjectMolecule.cpp:10004](layer2/ObjectMolecule.cpp#L10004), [CoordSet.h:145](layer2/CoordSet.h#L145) | `api/mmpymolx.py`, `jira/PYMOL-317.py`, `jira/PYMOL-1276.py` | Can be bridged with a Python callback (RDKit) after P-02. `text_type` needs an atom typer (e.g. MMFF types via RDKit) | M |
| C-04 | [ ] | `auto_copy_images`: copy to the clipboard automatically after `ray` / `draw` | [SceneRay.cpp:758](layer1/SceneRay.cpp#L758), [Scene.cpp:4083](layer1/Scene.cpp#L4083) | none | `cmd._copy_image` already works in the Qt GUI ([pymol_qt_gui.py:1170](modules/pmg_qt/pymol_qt_gui.py#L1170)). Only the C++ calls are missing: replace the `#ifdef` with `PParse(G, "cmd._copy_image(quiet=0)")` | S |
| C-05 | [ ] | CIF: `cif_metalc_as_zero_order_bonds` | [CifMoleculeReader.cpp:1948](layer2/CifMoleculeReader.cpp#L1948) | none | Small: when the setting is on, read `metalc` connections as zero-order bonds | S |
| C-06 | [ ] | `CmdM2ioFirstBlockProperties` (MAE first-block properties) | [Cmd.cpp:6337](layer4/Cmd.cpp#L6337) | `properties/*` | Part of F-01 | — |
| C-07 | [ ] | Shader antialiasing `antialias_shader` 1/2 (FXAA/SMAA) | `SettingInfo.h:796` | `settings/settings.py::testAA` | Run the test first; this may already work in open source, so it may only need un-gating | S–M |

## D. GUI items marked as Incentive-only

| ID | Status | Item | Location | Approach | Effort |
|---|---|---|---|---|---|
| G-01 | [ ] | Properties dialog: object-state user properties pane | [properties_dialog.py:80](modules/pmg_qt/properties_dialog.py#L80) | `get_property_list` / `set_property` already exist in [properties.py](modules/pymol/properties.py); add the tree category and editors | S |
| G-02 | [ ] | Properties dialog: atom-level user properties pane | [properties_dialog.py:86](modules/pmg_qt/properties_dialog.py#L86) | Same, using `set_atom_property` and `iterate p.*` | S |

## E. Behaviour differences noted in tests

| ID | Status | Difference | Source | Action |
|---|---|---|---|---|
| B-01 | [ ] | `fab ... dir=-1`: Incentive needs `sort` afterwards; open source is already sorted | [testing/tests/api/editor.py:50](testing/tests/api/editor.py#L50) | Probably fine as it is. Record it and move on |
| B-02 | [ ] | Version-dependent test expectations (`incentive: 1.8.4, open-source: 2.1`) | [testing/tests/api/querying.py:641](testing/tests/api/querying.py#L641) | Check that the pymolx version string doesn't break `requires_version` gating |
| B-03 | [ ] | `bg_image_filename` embedded as a `data:` URL when saving sessions | [testing/tests/jira/PYMOL-1571.py](testing/tests/jira/PYMOL-1571.py) | Embed the image on session save, decode it on load. Needs a GUI test run |
| B-04 | [ ] | `undo` after `remove` on discrete objects | [testing/tests/jira/PYMOL-1697.py](testing/tests/jira/PYMOL-1697.py) | Run the test with `--with-undo` first; it may already pass |
| B-05 | [x] | Multi-level undo/redo with Incentive's API (`undo_enable`, `undo_disable`, `undo`, `redo`, capability `multi_undo`). Session snapshots; one step per top-level command, per event-loop turn in the GUI; camera kept on undo. On by default in the GUI | [modules/pymolx/undo.py](modules/pymolx/undo.py) | The upstream undo tests ([testing/tests/undo/](testing/tests/undo/)) now run: 26 of 28 pass, the 2 others need `clean` (P-01) |

---

## F. Look and feel

These items are not marked in the code. We don't have access to
Incentive PyMOL, so the references are two annotated screenshots of the
Incentive PyMOL 3 main window from Schrödinger's documentation (dark
theme): a ligand-binding-site scene, and a session with the Scenes panel
and Timeline open. They are not committed here, as the images are
Schrödinger's. What they show:

- **Toolbar** under the menu: selection mode dropdown (pointer icon,
  "Residues", accent color), undo, redo, Zoom dropdown, Orient, Rock,
  Presets...; right-aligned: Builder..., Scenes, Draw/Ray dropdown
  (camera icon), `...` overflow.
- **Menu bar**: File, Edit, Build, Movie, Display, Setting, Scenes, Mouse,
  Wizard, Plugin, Help.
- **Content panel** (right of the viewer, called so in the docs): green
  dot for enabled entries, disabled ones greyed with near-invisible
  buttons, selections shown as `(name)`, rounded A/S/H/L/C buttons; the
  `All` row has purple A/S/H/L and a green C.
- **Toggle toolbar** under the content panel: Mouse (mode dropdown,
  "3-Button Viewing"), Wizard (wand), Sequence (`SEQ`), Timeline, Command
  (`>_`).
- **Command line** at the very bottom, full width, prompt `PyMOL >`; the
  output pane above it is toggled by `>_`; echoed commands in the accent
  color. There is no command line inside the viewer.
- **Scenes panel** (left dock, opened by Scenes): scene thumbnails named
  001, 002, ... each with a `...` menu; Save Scene and Add to Timeline
  buttons.
- **Timeline** (bottom dock, above the command line): composition tabs,
  ADD TRACK, a Camera track with keyframes, transport controls, Loop,
  time display, EDIT..., EXPORT COMP, SHOW INSPECTOR; a Timeline
  Inspector with camera position and rotation.

Anything not visible in the screenshot (menu contents, dialogs, light
theme) is our choice and marked *unconfirmed*. More reference
screenshots from Schrödinger's public docs would help most for: an open
A/S/H/L/C menu, the selection mode and Draw/Ray dropdowns, the `...` and
wand menus.

| ID | Status | Item | Where to work | Effort |
|---|---|---|---|---|
| L-01 | [x] | Application theme: dark Fusion palette plus stylesheet, purple accent. Colors, sizes and dividers measured from a running Incentive PyMOL 3.1.8 (neutral greys `#343434`/`#222222`, `#494949` dividers, `#9b81fd` accent). Like Incentive, the viewer shows the startup text until something is loaded or clicked (no splash image) | [modules/pymolx/gui/theme.py](modules/pymolx/gui/theme.py), [data/pymolx/styles/dark.qss](data/pymolx/styles/dark.qss) | M |
| L-02 | [~] | Our own icon set for the toolbar and menus. Done: pointer, undo, redo, camera, chevron | [data/pymolx/icons/](data/pymolx/icons/) | M |
| L-03 | [~] | Dock layout. Done: output pane and command line moved to the bottom, full width. Open: sequence viewer, scenes and movie panels; "Reset layout" action | [pymol_qt_gui.py](modules/pmg_qt/pymol_qt_gui.py), [scene_bin_gui.py](modules/pmg_qt/scene_bin_gui.py) | M |
| L-04 | [~] | Default settings at startup (shaders, AA, ray, cartoon, background), shipped as a startup module, not as `Setting.cpp` edits. **Mechanism done; `DEFAULTS` stays empty until values are confirmed side by side** | [modules/pymolx/defaults.py](modules/pymolx/defaults.py) | S |
| L-05 | [~] | Menu structure and wording matching Incentive. Done: `Scene` renamed to `Scenes` (top-level names now match the screenshot). Open: menu contents (*unconfirmed*) | [modules/pymol/_gui.py](modules/pymol/_gui.py) | M |
| L-06 | [ ] | Unified Preferences dialog to replace the raw list in [advanced_settings_gui.py](modules/pmg_qt/advanced_settings_gui.py) | new Qt dialog | M |
| L-07 | [ ] | Sequence viewer polish (Qt widget instead of the OpenGL strip) | new Qt widget | L |
| L-08 | [ ] | Welcome / start screen and recent-files list | new Qt widget | S |
| L-09 | [ ] | Builder: fragment library and editing workflow parity | [builder.py](modules/pmg_qt/builder.py) | M |
| L-10 | [ ] | Movie maker / timeline editing parity | [modules/pymol/movie.py](modules/pymol/movie.py), Qt panel | M |
| L-11 | [ ] | Mutagenesis wizard: rotamer library and preview parity | `modules/pymol/wizard/mutagenesis.py` | M |
| L-12 | [ ] | Go through the recent Incentive release notes and add any feature not listed here | this file | S |
| L-13 | [x] | Toolbar replacing the upstream grid of quick buttons. *Unconfirmed:* Presets... is a dropdown applying to all objects (the ellipsis may mean a dialog in Incentive); Zoom menu (All, Selection, Center, Reset View) and `...` menu (Get View, Unpick, Deselect, Properties, Rebuild, Movie controls) contents are ours | [modules/pymolx/gui/toolbar.py](modules/pymolx/gui/toolbar.py) | M |
| L-14 | [x] | Toggle toolbar: mouse mode dropdown (`pymol.menu.mouse_config`), Wizard menu, `SEQ` (`seq_view`), Timeline (open source's `movie_panel` until L-17 exists), `>_` output pane. The viewer's own command line is off (`internal_feedback=0`); typing in the viewer goes to the command line | [content_panel.py](modules/pymolx/gui/content_panel.py) | S |
| L-15 | [~] | Content panel replacing the OpenGL object panel (`internal_gui=0`). Done: rows from the core's panel list (new `_cmd.get_panel_list`, [panel.py](modules/pymolx/panel.py)), enabled dot, click to enable/disable, right-click for actions, A/S/H/L/C menus dispatched per object type as in the OpenGL panel, groups with expand/collapse, long names shortened with "...", menu labels in PyMOL's colors, a wizard panel (controls of the active wizard, e.g. mutagenesis) and a state bar (state stepping, rotamer strain). Open: drag to reorder, rename, multi-select | [content_panel.py](modules/pymolx/gui/content_panel.py), [layer3/Executive.cpp](layer3/Executive.cpp) | L |
| L-16 | [ ] | Scenes panel: left dock with scene thumbnails, a `...` menu per scene, Save Scene, Add to Timeline (upstream has a Qt scene panel in [scene_bin_gui.py](modules/pmg_qt/scene_bin_gui.py) to build on) | new Qt dock | M |
| L-18 | [x] | Superposition/Alignment dialog (Plugin menu): many-to-one (`extra_fit`) or one-to-one, method (align, super, cealign, usalign, fit), selections and states, alignment object, outlier rejection, command preview. Bundled plugin | [data/startup/alignment_gui/](data/startup/alignment_gui/) | S |
| L-19 | [x] | Interface Analysis (PISA-like), beyond Incentive: `interface_analysis` command and Plugin-menu dialog. Per chain pair: interface area, ΔiG (solvation parameters fitted to PDBe PISA: 0.66 kcal/mol mean error on held-out entries), H-bonds, salt bridges, disulfides, residue table, viewer display, CSV export | [modules/pymolx/interfaces.py](modules/pymolx/interfaces.py), [data/startup/interface_gui/](data/startup/interface_gui/), [tools/pymolx/calibrate_interfaces.py](tools/pymolx/calibrate_interfaces.py) | M |
| L-20 | [x] | Structure preparation and minimization, beyond Incentive: `fix_structure` (PDBFixer: missing atoms/residues, nonstandard residues, heterogens, hydrogens; new object, B-factors kept) and `minimize` (OpenMM: amber14/amber99sbildn/charmm36, GBn2/OBC2/vacuum, backbone or heavy-atom restraints, unselected atoms fixed, CUDA/OpenCL/CPU, one undo step); Structure Preparation dialog with background minimization and cancel. Optional deps: `pip install .[mm]` | [modules/pymolx/mm.py](modules/pymolx/mm.py), [data/startup/structure_prep_gui/](data/startup/structure_prep_gui/) | M |
| L-17 | [ ] | Timeline: compositions, tracks (camera, objects), keyframes, transport, export, inspector. Large; open source only has the OpenGL movie panel | new Qt dock | L |

## G. Bundled dependencies and packaging

| ID | Status | Item | Approach | Effort |
|---|---|---|---|---|
| D-01 | [ ] | APBS + pdb2pqr bundled, so the existing `apbs_gui` plugin works out of the box | conda packages in the installer | S |
| D-02 | [ ] | ffmpeg bundled for MP4 movie export | `imageio-ffmpeg`; point movie export at its binary | S |
| D-03 | [x] | Scientific stack: numpy, scipy, RDKit, gemmi, biopython | `pip install .[science]` ([pyproject.toml](pyproject.toml)); the `pymolx` command reports what is found | S |
| D-04 | [ ] | One-click installers for Linux, macOS and Windows with an app icon and file associations (`.pse`, `.pdb`, `.cif`, `.mae`) | conda `constructor` or briefcase | M |
| D-05 | [ ] | Pre-installed plugin set (lighting, APBS, …) checked on first launch | `modules/pmg_tk/startup/`, plugin manager | S |

## H. Branding and identity (required for the fork)

| ID | Status | Item | Location |
|---|---|---|---|
| H-01 | [x] | Our own splash image (pymolx), replacing `splash.png` / `ipymol.png`. Generated by [tools/pymolx/make_splash.py](tools/pymolx/make_splash.py); replace with final artwork when available | [commanding.py:325](modules/pymol/commanding.py#L325), `data/pymol/` |
| H-02 | [x] | About dialog: name, licence, contact (currently points to sales@schrodinger.com) | [pymol_qt_gui.py:905](modules/pmg_qt/pymol_qt_gui.py#L905) |
| H-03 | [x] | Help menu links to schrodinger.com pages: keep them as references or point them at pymolx docs | [_gui.py:869](modules/pymol/_gui.py#L869) |
| H-04 | [x] | Desktop app id `com.schrodinger.pymol` → our own id | [pymol_qt_gui.py:1218](modules/pmg_qt/pymol_qt_gui.py#L1218) |
| H-05 | [x] | Version string / product name in the startup banner | [modules/pymol/\_\_init\_\_.py:164](modules/pymol/__init__.py#L164), [layer1/Ortho.cpp](layer1/Ortho.cpp) |
| H-06 | [ ] | Package metadata: distribution name is still `pymol`, URLs and description point upstream. Renaming the distribution changes how existing installs upgrade, so decide deliberately | [pyproject.toml](pyproject.toml) |
| H-07 | [ ] | Legacy Tk GUI (`pmg_tk`) still shows upstream product names | [modules/pmg_tk/skins/normal/\_\_init\_\_.py:1279](modules/pmg_tk/skins/normal/__init__.py#L1279) |

All product naming comes from [modules/pymolx/branding.py](modules/pymolx/branding.py). The C startup banner duplicates it and must be kept in sync.

## I. Test infrastructure (do this first)

| ID | Status | Item | Approach | Effort |
|---|---|---|---|---|
| T-01 | [x] | Run the upstream Incentive tests for features we implement | An implementation calls `pymolx.features.register('P-04')`; [testing/pymolx_parity.py](testing/pymolx_parity.py) maps each parity ID to the upstream tests that verify it, and [testing/testing.py](testing/testing.py) runs a gated test once all its IDs are registered. Upstream test files stay unchanged. [tests/pymolx/test_foundation.py](testing/tests/pymolx/test_foundation.py) fails if a gated test is not mapped | S |
| T-02 | [ ] | Image-diff parity suite: render the same `.pml` scenes in Incentive and pymolx, compare with SSIM, track the score per scene in CI | `testing/parity/` (new) | M |
| T-03 | [ ] | Baseline: record which `requires('incentive')` tests pass after each milestone | CI job | S |

---

## Suggested order

1. **Foundation:** T-01, L-04, H-01…H-05, D-03.
2. **Quick wins:** C-04, P-04, P-06, P-07 (linear), F-03, F-04, F-05, C-05, G-01, G-02, P-02.
3. **High-value science:** P-01, P-03, P-09, C-03, F-01, F-02.
4. **Look and feel:** L-01…L-03, L-05, L-06, L-08, plus D-01, D-02, D-04.
5. **Core rendering:** C-01, C-02, P-07 (rigid-body morph), L-07.
