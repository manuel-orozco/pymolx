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
- Optional deps installed: rdkit, gemmi, scipy, openmm (CUDA), pdbfixer.

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
  `structure_prep_gui` (worker QThread pattern).
- Theme: colors only in `COLORS` (`modules/pymolx/gui/theme.py`), used as
  `$name` in `data/pymolx/styles/dark.qss` (literal `$` must be `$$`).
  Icons: `data/pymolx/icons/*.svg` (own artwork; `<name>-on.svg` = checked).
  Console font: `theme.CONSOLE_FONTS` (Consolas first, with fallbacks).
- Menus from `pymol.menu` data: `pymolx.gui.menus.PyMenu` + `fill_menu`
  (draws PyMOL color codes like `\900`).
- Main window parts: toolbar (`gui/toolbar.py`), content panel with
  wizard panel, state bar and toggle toolbar (`gui/content_panel.py`),
  mouse mode menu (`gui/mouse_modes.py`), sequence viewer dock above the
  viewer (`gui/sequence_viewer.py`, data in `pymolx/sequence.py`).
  Output/command line is the bottom dock `window.ext_window`.

## Gotchas found the hard way

- PySide6: `QMenu.addMenu("title")` submenus get deleted when the Python
  reference goes; create `QMenu(title, parent)` and `addMenu(menu)`.
  Fetching menus through temporary `action.menu()` wrappers can also fail
  in scripts; find live menus via `QApplication.allWidgets()`.
- `cmd.do` only queues commands (run later by the C loop).
- `iterate` variable for ATOM/HETATM is `type` (no `hetatm`).
- `fixed` is a selection keyword: never name an object "fixed".
- `get_area` on a selection counts the rest of the object; measure
  isolated copies.
- Settings changed during window setup become stored defaults
  (`reinitialize store` runs later in `adapt_to_hardware`).
- OpenMM picks `label_asym_id` as chain when it has more values than
  `auth_asym_id`; `pymolx.mm` writes its own mmCIF to avoid this.
- `_cmd.get_seq_change_counts` (changed, dirty) tells when residues,
  colors or selections may have changed. Queries on a selection
  expression create temporary selections and bump "dirty" themselves:
  re-read the counters after your own queries.
- Kill test processes by PID, never `pkill -f` with a pattern that
  matches your own shell command.

## Working with the owner

- Explain in plain language; show what was verified (tests, screenshots).
- Don't commit or push unless asked; the owner usually commits to
  `master` directly. Keep build artifacts and downloads out of the repo.
- Decisions already made: keep ESC text/graphics toggle and the startup
  text in the viewer (like Incentive); Display > Background has Dark Navy
  `#000430` and Custom...; `interface_analysis` ΔiG is calibrated to PISA;
  4× multisampling by default (`options.multisample`, `-E 0` = off).
