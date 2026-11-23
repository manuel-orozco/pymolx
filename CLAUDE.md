# pymolx

A fork of Open-Source PyMOL that aims to look and act like Incentive
(commercial) PyMOL, plus extra tools. Owner: Manuel Orozco
(github.com/manuel-orozco/pymolx). Status of every parity item, with IDs
used in code and commits: [INCENTIVE_PARITY.md](INCENTIVE_PARITY.md).

## Build and run

- Python env: `../pymolx-env` (not activated; call binaries by path).
- Build + install after **any** change (Python, data or C++; incremental,
  ~15 s). Without `use-msgpackc=no` the build fails on a missing header:

      ../pymolx-env/bin/pip install -v --no-build-isolation --no-deps \
          --config-settings use-msgpackc=no . > /tmp/build.log 2>&1

- Run: `../pymolx-env/bin/pymol` (GUI). Headless: `pymol -cq script.py`.
- Optional deps installed: rdkit, gemmi, scipy, openmm (CUDA), pdbfixer,
  propka, freesasa (aggregation prediction), anarci + hmmer (antibody
  numbering), torch (CUDA),
  parapred-pytorch and Paragraph (paratope prediction, both installed
  from GitHub; parapred's wheel leaves out its weights, which were
  fetched into its package's weights/ folder by hand).
- torch is 2.15.0+cu132, the plain PyPI wheel (`pip install
  torch==<version>`; the pytorch.org CUDA indexes lag a release or
  two). An installed `+cpu` wheel satisfies `torch==<version>`, so pip
  skips the swap: uninstall torch first. The machine has two GPUs,
  device 0 the 4090 and device 1 the 5090, and `torch.device("cuda")`
  (what ParaSurf asks for) takes device 0.
- ParaSurf (the third paratope predictor) is not a package: it is
  cloned at `/media/pac/data/apps/parasurf`, beside the other tools
  (`paratope.parasurf_home` looks for "ParaSurf" or "parasurf" next to
  the environment, next to its checkout and in the home directory, or
  takes `$PARASURF_HOME`), with its three Paragraph-expanded `.pth`
  weight files downloaded into its `weights/`. It also needs
  openbabel-wheel, torchsummary, scikit-learn, jsonpickle, matplotlib
  and h5py (installed), and the DMS surface program, built from its
  `dms/` into the environment with
  `make LIBDIR=$ENV/lib BINDIR=$ENV/bin MANDIR=$ENV/share/man/man1
  install`, with ENV=/media/pac/data/apps/pymolx/pymolx-env (no sudo;
  the shipped binary wants /usr/local/lib/dms, and `parasurf_dms` also
  looks in the environment's bin). An Fv of a Fab (1n8z) takes 16 s
  at surface detail 0.1 and 33 s at ParaSurf's own 0.3, against 86 s
  on the CPU, for the same scores to three decimals; a run is still
  not instant, because ANARCI, DMS and pdb2pqr come first.

- Aggrescan4D is re-implemented in `pymolx.aggrescan` from the
  published Aggrescan3D code, cloned at `/media/pac/data/apps/aggrescan`
  (`aggrescan3d/`, MIT; see its README.md for what replaces what). The
  score, its constants and the pH scale are ported into the module, so
  the clone is reference only. What A4D does with paid tools is done
  with OpenMM + PDBFixer (optimisation), EvoEF (mutant dG), implicit
  solvent dynamics (dynamic mode, in place of CABS-flex), PROPKA (3D
  pKa, in place of pKa-ANI) and a pLDDT cutoff (in place of
  AlphaCutter).
- The port is checked against real Aggrescan3D output in
  `/media/pac/data/apps/aggrescan3d/outputs` (`test_aggrescan.py`
  skips if absent). Matching it needs the relative exposure rounded to
  one decimal (`RSA_DECIMALS`), which is the precision A3D reads back
  from freesasa's `.rsa` text: a residue at 9.96% rounds to 10.0 and
  crosses the buried threshold, which moves its score and its
  neighbours'. Without the rounding a handful of residues per
  structure differ by up to 0.6.
- Do not score a mutation by the difference of two OpenMM potential
  energies: the mutant and the wild type differ in composition and
  that swamps any stability effect (Leu to Arg came out at -183
  kcal/mol). Minimizing with everything but the new side chain frozen
  is no better (thousands of kcal/mol: the rebuilt side chain cannot
  escape a clash). EvoEF's `ComputeStability` against the `_WT`
  companion it writes beside each mutant gives a usable dG; it takes
  FoldX's own mutant-file syntax, which A3D also uses.
- The "aggregation" preset lives in `aggrescan.preset`, hooked into
  upstream from `preset.aggregation` (`modules/pymol/preset.py`) and
  listed in `menu.presets` (`modules/pymol/menu.py`). The toolbar's
  Presets button builds from that same `pymol.menu` data on
  `aboutToShow`, so a test has to emit that signal before the menu has
  any actions (showing it would block). Set a preset's transparency on
  the object, not on the selection `_prepare` hands back: an atom level
  transparency quietly overrides a later
  `set transparency, ..., obj`.
- Scores are written into a PDB the way A3D writes output.pdb: the
  B-factor column (61-66) at `%6.2f`, occupancy untouched. Pad the
  line to 80 characters before slicing it, or the element symbol at
  the end of a PyMOL-written line is lost.
- A 5 A sphere does not simply lower the scores: it raises them,
  because most neighbours within 10 A are polar and pull downwards.
  What it does is keep a residue's score near its own propensity.

## Tests

- pymolx tests (pytest, run inside PyMOL):

      cd testing && ../../pymolx-env/bin/pymol -cq -d 'python
      import sys, pytest, pymol
      pymol.__path__.append(".")
      pytest.main(["-q", "--tb=short", "-p", "no:cacheprovider", "tests/pymolx"])
      sys.stdout.flush()
      python end'

- Full suite: `cd testing && ../../pymolx-env/bin/pymol -cq runall.pml`
  (unittest phase, then pytest phase). Known failures, all environment
  or not-yet-implemented: unittest 5 errors (MMTF x4, glTF), pytest 8
  (bcif x4, protonate_fallback x2 = cwd issue, test_undo_clean x2 = needs
  `clean`, P-01). Anything else failing is a regression.
- Upstream tests marked `requires('incentive')` run once their parity ID
  is registered (`pymolx.features.register('P-04')`, map in
  `testing/pymolx_parity.py`).
- GUI checks: script with `QtCore.QTimer.singleShot` steps, find the main
  window via `QApplication.topLevelWidgets()`, `widget.grab().save(...)`
  and look at the PNG. Don't grab the whole screen.

## Layout and conventions

- New code goes in `modules/pymolx/` (Qt-free core) and
  `modules/pymolx/gui/` (Qt). Upstream files get minimal hooks, commented
  `# pymolx: ...`. Match surrounding style: `_self=cmd` keyword, PyMOL
  docstring format (DESCRIPTION / USAGE / ARGUMENTS) for commands.
- Commands: register in `pymolx._init` (`modules/pymolx/__init__.py`)
  with `_self.extend(name, func)` **and** `setattr(_self, name, func)`.
  State-changing commands also go in `UNDOABLE` (`modules/pymolx/undo.py`)
  and must be registered before `undo.install`.
- Plugins: `data/startup/<name>/__init__.py` with `__init_plugin__` +
  `plugins.addmenuitemqt(...)`; dialog in `main.py`; reach the window with
  `pymol.gui.get_qtwindow()`. Examples: `alignment_gui`, `interface_gui`,
  `structure_prep_gui` (worker QThread pattern), `cluster_gui` (RMSD
  clustering of e.g. antibody-antigen models, `pymolx.cluster`),
  `compare_gui` (similarity coloring of close relatives, `pymolx.compare`),
  `seqcompare_gui` (the same by sequence, `pymolx.seqcompare`),
  `pocket_gui` (residues near a ligand, `pymolx.pocket`),
  `antibody_gui` (full IgG from Fab/scFv/VHH arms, hinge, Fc and
  C-terminal fusions, `pymolx.antibody`),
  `prediction_gui` (FASTA plus modifications, glycans and ligands as an
  AlphaFold 3 / RoseTTAFold3 / ESMFold2 input file, `pymolx.prediction`),
  `epitope_gui` (the same for antibody-antigen cofolding with AlphaFold 3:
  chain roles, many seeds, samples per seed, the antibody's MSA;
  `pymolx.epitope`, subclasses `PredictionDialog` and replaces its
  `_settings_box`),
  `cdr_gui` (CDR loops of an antibody under IMGT/Kabat/Chothia/Martin/AHo
  numbering, from a structure or a sequence, `pymolx.cdr`),
  `paratope_gui` (binding residues by Parapred, Paragraph or ParaSurf,
  `pymolx.paratope`; ParaSurf runs in a worker QThread and PyMOL is
  only called before and after it),
  `aggregation_gui` (aggregation propensity and its pH dependence,
  `pymolx.aggrescan`; the whole analysis is a `Job` the dialog runs in
  a worker QThread and polls for `job.status`, as `structure_prep_gui`
  does; `Job.save` writes the run out as CSV and as scored PDB files,
  one per pH, so the job's temporary folder is kept until the next run
  or the dialog closes).
- Theme: colors only in `COLORS` (`modules/pymolx/gui/theme.py`), used as
  `$name` in `data/pymolx/styles/dark.qss` (literal `$` must be `$$`).
  Icons: `data/pymolx/icons/*.svg` (own artwork; `<name>-on.svg` = checked).
  Font size: `theme.FONT_SIZE` (11 pt) for all GUI text: application
  font, `$font_size` in the stylesheet, console font. Console font
  (output, command line, sequence viewer): `theme.console_font()`,
  first of `theme.CONSOLE_FONTS` (DejaVu Sans Mono), sized so its letters
  are as tall as the interface font's. No bold text except titles (and
  the user's choice of "Courier New Bold" in Display > Sequence > Font).
- Menus from `pymol.menu` data: `pymolx.gui.menus.PyMenu` + `fill_menu`
  (draws PyMOL color codes like `\900`).
- Main window parts: toolbar (`gui/toolbar.py`), content panel with
  wizard panel, state bar and toggle toolbar (`gui/content_panel.py`),
  mouse mode menu (`gui/mouse_modes.py`), sequence viewer dock above the
  viewer (`gui/sequence_viewer.py`, data in `pymolx/sequence.py`).
  Output/command line is the bottom dock `window.ext_window`. "Scripts"
  in the menu bar (`gui/scripts.py`) is File > Run Script... starting in
  `scripts/` of the source
  checkout: `setup.py` writes its path to the installed
  `pymolx/_source.py` (`pymolx.paths`).

## Gotchas found the hard way

- PySide6: `QMenu.addMenu("title")` submenus get deleted when the Python
  reference goes; create `QMenu(title, parent)` and `addMenu(menu)`.
  Fetching menus through temporary `action.menu()` wrappers can also fail
  in scripts; find live menus via `QApplication.allWidgets()`.
- PySide6: `signal.connect` returns a Connection, so `@signal.connect`
  as a decorator replaces the function; connect after defining it.
- PySide6: an attribute named like a virtual Qt method (e.g.
  `self.metric`, QPaintDevice::metric) replaces it: Qt calls it while
  painting and crashes (segfault on `grab()`).
- `cmd.do` only queues commands (run later by the C loop).
- `iterate` variable for ATOM/HETATM is `type` (no `hetatm`).
- Many words are selection keywords (`fixed`, `gap`, `none`, `alt`, `b`,
  `first`, `bonded`, ...): never name objects so (PyMOL appends "_").
- `iterate` prints the value of an expression (e.g. `d.setdefault(...)`):
  write an assignment (`_ = ...`); and `index` is an atom property there,
  don't use it as the space dict's name.
- `get_area` on a selection counts the rest of the object; measure
  isolated copies.
- Settings changed during window setup become stored defaults
  (`reinitialize store` runs later in `adapt_to_hardware`).
- OpenMM picks `label_asym_id` as chain when it has more values than
  `auth_asym_id`; `pymolx.mm` writes its own mmCIF to avoid this.
  `pymolx.prep` gives waters/ligands chains of their own (else sequences,
  caps and terminal H go wrong) and maps them back on loading.
- OpenMM topologies keep bonds of standard residues only: ligand bonds
  are copied back from PyMOL (`mm._het_bonds`, `mm._restore_bonds`).
  `ForceField.getUnmatchedResidues` ignores CHARMM terminal patches:
  only leave out non-polymer residues.
- OpenMM minimization (handled in `prep.PrepJob._minimize`): on CUDA,
  nothing moves after `context.reinitialize()` (use a new context); the
  tolerance is an RMS over all forces, so big systems stop with strained
  atoms (restart while the largest force is high). PDBFixer builds
  missing residues as D-amino acids or cis peptides at times (mirrored,
  then guided by torsion restraints); `1 - cos` guides don't push an
  exactly cis peptide (wrapped harmonic ones do).
- `_cmd.get_seq_change_counts` (changed, dirty) tells when residues,
  colors or selections may have changed. Queries on a selection
  expression create temporary selections and bump "dirty" themselves:
  re-read the counters after your own queries.
- `sendEvent` mouse events don't drag `QMainWindow` separators in the
  real GUI (they do offscreen); drive such drags with `QTest`.
- Plugin preferences (`~/.pymolpluginsrc.py`, `plugins.pref_get`) are
  read after the main window is built: read them again from a
  `QTimer.singleShot(0, ...)` set up during window setup.
- `cmd.load` is wrapped twice: `pymolx.missing` (adds PDB REMARK 465
  residues as CA atoms without coordinates, like the mmCIF reader) inside
  `pymolx.undo` (one undo step). Wrap commands before `undo.install` and
  update `_self.keyword` entries too. Atoms without coordinates: selection
  keyword `present` excludes them; `count_atoms` includes them.
- `load` also runs scripts (`.py`, `.pml`): code wrapping it must not
  take a script for the source of the objects it creates
  (`pymolx.files`).
- Installing while the GUI runs leaves it half old: `pymolx._init`
  imports the core modules at startup, so they stay as they were, while
  a plugin dialog is imported the first time it is opened and comes
  from the new file ("module 'pymolx.paratope' has no attribute
  'PARASURF_MESH'"). Restart PyMOL after a build.
- Never `import pymolx...` inside `execapp()` (pymol_qt_gui): it makes
  `pymolx` a local name there and the GUI doesn't start; import at the
  top of the file.
- Qt delivers mouse moves without a pressed button only to widgets with
  `setMouseTracking(True)`; `QTest.mouseMove` may not reach the real GL
  view, send `QMouseEvent`s to it instead.
- GUI test scripts run with the owner's real plugin preferences: swap in
  `MemoryPrefs()` before toggling anything that saves a preference.
- Picking (`_cmd.pick_at`) draws atoms in code colors (mostly red) into
  the frame buffer: run it only in the GL widget's `paintGL`
  (`pymolwidget.before_draw` hooks), else Qt may show it when another
  widget repaints (red flashes). Inside `with widget:` all the same:
  PyMOL draws nothing without its valid-context flag, even in paintGL.
- A new mouse button action (`layer1/ButMode.h`, `but_act_code` in
  `modules/pymol/controlling.py`) has to be handled in three places in
  `layer1/SceneMouse.cpp`: the switch that fills `selName`
  (`ExecutiveGetActiveSeleName`), the switch that acts on the click, and
  the `sel_mode_kw` branch of `SceneClick`. Miss one and the click binds
  but does nothing. In a GUI script, a click's effect shows up a step
  later (the pick happens in the next draw), and `cmd.mouse` into an
  editing mode deselects.
- PyMOL orders a chain's atoms by `segi` first: different segments
  within a chain scramble its residue order (and exports of it). Tag
  residues with the `custom` atom property instead (selectable:
  `custom HNGE`, kept in sessions).
- `alter ..., resv += n` drops insertion codes (28A becomes a second
  28); renumber with `resi = str(resv + n) + resi.lstrip("-0123456789")`.
- Objects named like selection keywords (`model`, `chain`, ...) get a
  `_` appended by PyMOL; the check isn't exposed to Python
  (`antibody.object_name` probes it with a one-atom object).
- Qt reads an `<a href="A:12">` as a URL whose *scheme* is "A", and
  lower-cases it: anchors carrying data need a fragment (`href="#A:12"`,
  read back with `url.fragment()`).
- A dialog a test only closes with `deleteLater()` is still alive when
  PyMOL quits (nothing runs the queued deletion in `pymol -cq`), and the
  C layer's teardown then segfaults, after every test has passed. Destroy
  test dialogs at once: `shiboken6.delete(dialog)`. The segfault shows up
  only in the full suite (`PYTHONFAULTHANDLER=1` points at `internal._quit`),
  not in the plugin's own test file, and never in the real GUI.
- A `QMenu.exec` in a slot blocks: a test that emits the signal hangs.
  Build the menu in a method the test can call, and disconnect the
  dialog's own slot before emitting.
- A full-IgG template's Fab arms lie where the Fc's partners bind, so
  a receptor carried in with the Fc (`fc_antigens`) lands inside an arm
  whatever template is used; 1HZH, 1IGT and 1IGY all put more of an arm
  there than a refined model does, and no intact-IgG structure has its
  arms posed for a bound receptor (3AY4 and the like are Fc + receptor
  only, the Fabs being flexible and unresolved). `swing_arms=1` turns
  the arm at its hinge instead. Turn it about the backbone *carbon*
  that joins the arm to the hinge, not the alpha carbon: that keeps the
  peptide bond's length exactly, where pivoting on CA stretches it
  (1.33 -> 1.68 A) and needs a relaxation afterwards.
- The antibody builder numbers its stitched chains 1..N, but chains
  carried in with a part (antigens, receptors) keep their own numbering:
  use `assembly.numbers(chain)`, not `enumerate(residues, 1)`.
- Renumbering residues one at a time (`alter ... resi = N` per residue)
  merges them when a new number collides with one not yet renumbered
  (e.g. a chain numbered 1..470 renumbered to IMGT 1..128). Build a
  {(chain, resi): (chain, resi)} map and do one pass:
  `alter obj, "chain, resi = NEW.get((chain, resi), (chain, resi))"`.
- ANARCI needs the real HMMER binaries, not the PyPI `hmmer` wheel's
  console shims: the shim's temporary name reaches hmmscan's output
  header, where Biopython's parser expects "hmmscan", so every alignment
  coordinate comes back None and ANARCI dies in `_domains_are_same`
  ("'<' not supported between instances of 'NoneType'"). Pass
  `hmmerpath=cdr.hmmer_path()`. Not a Biopython version problem (1.85
  fails the same way).
- Kill test processes by PID, never `pkill -f` with a pattern that
  matches your own shell command.
- Quitting the GUI by closing the window, File > Quit, the viewer
  menu's quit or typing `quit` asks to save a changed session (`pymolx.unsaved`); a script's
  `cmd.quit()` doesn't, so GUI test scripts still end with it.

## Working with the owner

- Explain in plain language; show what was verified (tests, screenshots).
- Don't commit or push unless asked; the owner usually commits to
  `master` directly. Keep build artifacts and downloads out of the repo.
- Decisions already made: keep the ESC text/graphics toggle, but no
  startup text in the viewer (`branding.SHOW_STARTUP_TEXT`); Display > Background has Dark Navy
  `#000430` and Custom...; `interface_analysis` ΔiG is calibrated to PISA;
  Display > Quality > Maximum Quality is the default (`defaults.QUALITY`);
  stick_radius 0.15 and line_width 1.0 px by default (`defaults.DEFAULTS`;
  thinner cylinder lines break up into dashes); selected atoms shown as
  light cyan dots (transparency 0.2) like Maestro (`selection_round_points` on; round
  indicator texture changed in `layer3/Executive.cpp`);
  no default multisampling (tried, didn't help); new objects get grey
  carbons (`grey70`, via `auto_color` 0 in `defaults.DEFAULTS` and
  `AtomInfoUpdateAutoColor`; the named color "carbon" stays green) and
  white hydrogens (`AtomInfoGetColor`, also with `auto_color` 0; "hydrogen"
  stays 0.9 grey), also the first entry of Color > by element (its
  "HNOS..." entry, `util.cnc`, removed); the output pane starts hidden
  (`content_panel.SHOW_OUTPUT`), only the command line shows until ">_";
  Structure Preparation is a Maestro-like workflow (`pymolx.prep`,
  `pymolx.hbonds` for flips/rotatable H, `pymolx.diagnostics`);
  restraints hold input atoms only; MolProbity's Top8000 grids
  (CC BY 4.0, credit in `data/pymolx/molprobity/README.md`) give
  Ramachandran/rotamer outliers matching wwPDB reports; metal sites
  are held fixed in the minimization (OpenMM's nonbonded ions drift away
  or, without van der Waals, collapse onto thiolates);
  Display > Lighting > Softer is the startup lighting (`lighting.DEFAULT_PRESET`
  in `defaults.DEFAULTS`; "lighting default" = softer, PyMOL's own is
  "lighting classic", menu "Classic PyMOL"); 3-Button Viewing binds
  middle drag Slab, right drag Move, wheel IMvZ = zoom in on scroll up
  (`controlling.mode_dict`;
  dragging Slab is pymolx's addition in `layer1/SceneMouse.cpp`, up =
  thicker).
