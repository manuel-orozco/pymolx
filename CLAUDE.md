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
  propka.

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
  C-terminal fusions, `pymolx.antibody`).
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
- PyMOL orders a chain's atoms by `segi` first: different segments
  within a chain scramble its residue order (and exports of it). Tag
  residues with the `custom` atom property instead (selectable:
  `custom HNGE`, kept in sessions).
- `alter ..., resv += n` drops insertion codes (28A becomes a second
  28); renumber with `resi = str(resv + n) + resi.lstrip("-0123456789")`.
- Objects named like selection keywords (`model`, `chain`, ...) get a
  `_` appended by PyMOL; the check isn't exposed to Python
  (`antibody.object_name` probes it with a one-atom object).
- Kill test processes by PID, never `pkill -f` with a pattern that
  matches your own shell command.

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
