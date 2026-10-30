'''
Antibody Builder dialog: a schematic antibody whose parts (two arms, each
a Fab, an scFv, a VHH or none; the hinge; the Fc) are clicked to give
their structures; a full-length template; the model built by
pymolx.antibody (placement in the GUI thread, refinement in a worker
thread).
'''

import traceback

from pymol import cmd
from pymol.Qt import QtCore, QtGui, QtWidgets

from pymolx import antibody, prep
from pymolx.gui import theme

Qt = QtCore.Qt

_dialog = None

# schematic, in logical coordinates (scaled to the widget)
SIZE = QtCore.QSizeF(440.0, 500.0)
CX = 220.0
ARM_ANGLE = 42.0           # degrees from vertical
ARM_BASE_Y = 200.0
BOX_W, BOX_H, BOX_GAP = 30.0, 62.0, 6.0
CHAIN_X = 18.0             # column offset from the arm/Fc axis

FORMAT_ITEMS = [('Fab', 'fab'), ('scFv', 'scfv'),
                ('VHH (single domain)', 'vhh'), ('None (one-armed)', 'none')]
ARM_CHAIN_LABELS = {'fab': ('Heavy chain', 'Light chain'),
                    'scfv': ('scFv chain',), 'vhh': ('VHH chain',),
                    'fusion': ('scFv or VHH chain',)}
FUSIONS = tuple('fusion_' + site for site in antibody.FUSION_SITES)
ORIENTATION_ITEMS = [('As in the structure', ''), ('VH-VL', 'vh-vl'),
                     ('VL-VH', 'vl-vh')]
SHAPE_ITEMS = [('Relaxed (near the chain)', 'relaxed'),
               ('Extended', 'extended'),
               ('Extended, then collapsed by dynamics', 'collapsed')]
POSITIONS = [('On the template\'s Fv', 'fv'),
             ('Next to the hinge (in place of CH1/CL)', 'hinge')]
ARM_DESCRIPTIONS = {
    'fab': 'A Fab: heavy chain (VH-CH1) and light chain (VL-CL), e.g. a '
           'crystal structure or a predicted model.',
    'scfv': 'An scFv: VH and VL in one chain (either order). It is placed '
            'next to the hinge, like the Fv of an scFv-Fc; tags after the '
            'last domain are left out.',
    'vhh': 'A single-domain antibody (VHH, nanobody), placed next to the '
           'hinge.',
    'fusion': 'An scFv or VHH appended to the C-terminus of the chain '
              'through a linker: placed pointing away from the rest, the '
              'linker built by the refinement.',
}
CHAIN_LABELS = {
    'fc': ('Chain 1', 'Chain 2'),
    'hinge': ('Chain (Fab1 side)', 'Chain (Fab2 side)'),
    'template': ('Heavy chain 1', 'Light chain 1', 'Heavy chain 2',
                 'Light chain 2'),
}
DESCRIPTIONS = {
    'fc': 'The Fc: two heavy-chain fragments (CH2-CH3), with or without '
          'the lower hinge; glycans are kept.',
    'hinge': 'The hinge joins each arm to an Fc chain. By '
             'default the template\'s hinge is used; give a sequence (e.g. '
             'another isotype) or a structure instead.',
    'template': 'A full-length IgG (two heavy, two light chains): the parts '
                'are superposed on it and its hinge joins them.',
}
TITLES = dict(antibody.PART_LABELS, template='Template',
              fusions='Fusions', **{
                  'fusion_' + site: 'Fusion on %s' % name
                  for site, name in antibody.FUSION_SITES.items()})
ARMS = ('fab1', 'fab2')


def _color(hex_color, lighter=0.0):
    r, g, b = antibody._lighter(hex_color, lighter)
    return QtGui.QColor.fromRgbF(r, g, b)


def _linker_text(linker):
    '''
    "GGGGSGGGGS" -> "(G4S)2"; "" -> "none"
    '''
    if not linker:
        return 'none'
    if len(linker) % 5 == 0 and linker == 'GGGGS' * (len(linker) // 5):
        return '(G4S)%d' % (len(linker) // 5)
    return linker


def summary(part, spec, fmt='fab', same_as=None):
    '''
    One line for a part's structure (spec: PartDialog.spec()); arms: of
    format fmt; same_as: the arm whose structure an empty one repeats
    '''
    if part.startswith('fusion_') and spec:
        return '%s %s (chain %s%s), linker %s' % (
            antibody.FORMAT_LABELS.get(spec.get('kind'), ''),
            spec['object'], spec['chains'][0],
            ', ' + spec['orientation'].upper()
            if spec.get('orientation') else '',
            _linker_text(spec['linker']) +
            ({'extended': ', extended',
              'collapsed': ', collapsed'}.get(spec.get('linker_shape'), '')))
    if part in ARMS and fmt == 'none':
        return 'none (one-armed)'
    if spec is None and same_as:
        return 'same as %s' % same_as
    if spec is None:
        return 'click to give its structure' if part != 'hinge' else \
            'from the template'
    if part == 'hinge':
        if spec['mode'] == 'template':
            return 'from the template'
        if spec['mode'] == 'sequence':
            text = spec['sequence']
            return 'sequence %s%s (%d aa)' % (
                text[:12], '...' if len(text) > 12 else '', len(text))
    chains = spec.get('chains')
    if part == 'template' and not chains:
        known = dict(antibody.TEMPLATES).get(spec['object'], '')
        return '%s%s' % (spec['object'],
                         ': ' + known.split(',')[0] if known else '')
    if part in ARMS and chains and fmt == 'fab':
        return '%s (heavy %s, light %s)' % (spec['object'], chains[0] or '-',
                                            chains[1] or '-')
    if part in ARMS and chains:
        return '%s (chain %s%s)%s' % (
            spec['object'], chains[0] or '-',
            ', ' + spec['orientation'].upper()
            if spec.get('orientation') else '',
            ', linker %s' % _linker_text(spec['linker'])
            if spec.get('linker') else '')
    if part == 'template' and chains:
        return '%s (%s/%s, %s/%s)' % ((spec['object'],) + tuple(chains))
    return '%s (%s)' % (spec['object'], ', '.join(c or '-' for c in chains))


class Schematic(QtWidgets.QWidget):
    '''
    The IgG drawn as boxes (domains) and lines (hinge); clicking a part
    emits clicked(part)
    '''
    clicked = QtCore.Signal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setObjectName('antibody_schematic')
        self.setMouseTracking(True)
        self.setMinimumSize(330, 360)
        self.setSizePolicy(QtWidgets.QSizePolicy.Policy.Expanding,
                           QtWidgets.QSizePolicy.Policy.Expanding)
        self.assigned = {}     # part: bool
        self.texts = {}        # part: summary
        self.hover = ''
        self.formats = {'fab1': 'fab', 'fab2': 'fab'}
        self.titles = dict(antibody.PART_LABELS)
        self.fusions = {}      # site part: unit kind (scfv, vhh)
        self.position = 'fv'   # of scFv and VHH arms
        self._build_shapes()

    # geometry

    def _box(self, transform, x, y):
        path = QtGui.QPainterPath()
        path.addRoundedRect(QtCore.QRectF(x - BOX_W / 2, y - BOX_H / 2,
                                          BOX_W, BOX_H), 6, 6)
        return transform.map(path)

    def _fusion_unit(self, site, transform, side, end):
        '''
        A fused unit (or a dashed site to add one, heavy chains only) at
        transform, joined to the chain end by a linker line
        '''
        kind = self.fusions.get(site)
        top = transform.map(QtCore.QPointF(0, -BOX_H / 2))
        if not kind:
            if site.startswith('fusion_heavy'):
                path = QtGui.QPainterPath()
                path.addRoundedRect(QtCore.QRectF(-BOX_W / 2, -BOX_H / 2,
                                                  BOX_W, BOX_H), 6, 6)
                self.sites.append((site, transform.map(path),
                                   transform.map(QtCore.QPointF(0, 0))))
            return
        if kind == 'scfv':
            for kind_x, chain_kind in ((-side * 16, 'heavy'),
                                       (side * 16, 'light')):
                self.boxes.append((site, chain_kind,
                                   self._box(transform, kind_x, 0)))
            top = transform.map(QtCore.QPointF(-side * 16, -BOX_H / 2))
        else:
            self.boxes.append((site, 'heavy', self._box(transform, 0, 0)))
        line = QtGui.QPainterPath(end)
        line.lineTo(top)
        self.linker_lines.append((site, line))

    def _build_shapes(self):
        '''
        self.boxes: [(part, chain kind, path)]; self.hinge_lines: [path];
        self.label_at: {part: point}; self.sites: [(fusion part, path,
        center)] where a fusion can be added; self.linker_lines
        '''
        self.boxes = []
        self.label_at = {}
        self.sites = []
        self.linker_lines = []
        constant = -(BOX_H / 2 + 12)
        variable = constant - BOX_H - BOX_GAP
        heavy_bases = []
        for part, side in (('fab1', -1), ('fab2', 1)):
            fmt = self.formats[part]
            t = QtGui.QTransform()
            t.translate(CX + side * 12, ARM_BASE_Y)
            t.rotate(side * ARM_ANGLE)
            # heavy chain (VH) inside (toward the axis), light chain (VL)
            # outside; an scFv or VHH where a Fab's Fv is (a linker line to
            # the hinge) or next to the hinge
            heavy_x, light_x = -side * CHAIN_X, side * CHAIN_X
            single = (variable,) if self.position == 'fv' else (constant,)
            rows = {'fab': (constant, variable), 'scfv': single,
                    'vhh': single, 'none': ()}[fmt]
            if fmt in ('scfv', 'vhh') and self.position == 'fv':
                line = QtGui.QPainterPath(t.map(QtCore.QPointF(heavy_x,
                                                               -12)))
                line.lineTo(t.map(QtCore.QPointF(heavy_x,
                                                 variable + BOX_H / 2)))
                self.linker_lines.append((part, line))
            for y in rows:
                self.boxes.append((part, 'heavy', self._box(t, heavy_x, y)))
                if fmt != 'vhh':
                    self.boxes.append((part, 'light',
                                       self._box(t, light_x, y)))
            heavy_bases.append(t.map(QtCore.QPointF(heavy_x, -12)))
            top = rows[-1] if rows else constant + BOX_H / 2
            self.label_at[part] = t.map(QtCore.QPointF(0, top - 70))
            # a unit fused to the light chain: beside its constant domain
            site = 'fusion_light%d' % (1 if side < 0 else 2)
            if fmt == 'fab' and self.fusions.get(site):
                f = QtGui.QTransform(t)
                f.translate(light_x + side * (BOX_W + 34), -34)
                self._fusion_unit(site, f, side,
                                  t.map(QtCore.QPointF(light_x, -12)))
        fc_top = []
        for side in (-1, 1):
            t = QtGui.QTransform()
            t.translate(CX + side * CHAIN_X, 300)
            t.rotate(-side * 6)
            ch2 = -(BOX_H + BOX_GAP) / 2
            ch3 = (BOX_H + BOX_GAP) / 2
            for y in (ch2, ch3):
                self.boxes.append(('fc', 'heavy', self._box(t, 0, y)))
            fc_top.append(t.map(QtCore.QPointF(0, ch2 - BOX_H / 2)))
            # a unit fused to the heavy chain (or a site to add one) below
            site = 'fusion_heavy%d' % (1 if side < 0 else 2)
            f = QtGui.QTransform()
            f.translate(CX + side * 44, 300 + ch3 + BOX_H + 22)
            f.rotate(side * 18)
            self._fusion_unit(site, f, side,
                              t.map(QtCore.QPointF(0, ch3 + BOX_H / 2)))
        # hinge: one line per heavy chain, two disulfide bars
        self.hinge_lines = []
        for k, side in enumerate((-1, 1)):
            path = QtGui.QPainterPath(heavy_bases[k])
            path.lineTo(CX + side * 8, ARM_BASE_Y + 14)
            path.lineTo(CX + side * 8, ARM_BASE_Y + 40)
            path.lineTo(fc_top[k])
            self.hinge_lines.append(path)
        for y in (ARM_BASE_Y + 22, ARM_BASE_Y + 31):
            path = QtGui.QPainterPath(QtCore.QPointF(CX - 8, y))
            path.lineTo(CX + 8, y)
            self.hinge_lines.append(path)
        self.hinge_area = QtCore.QRectF(CX - 34, ARM_BASE_Y - 4, 68, 56)
        # hinge and Fc labels left of them, or right (a unit fused to the
        # left light chain), or at the bottom left (fused to both)
        left = not self.fusions.get('fusion_light1')
        right = not self.fusions.get('fusion_light2')
        self.label_side = 'left' if left else 'right' if right else 'bottom'
        if self.label_side == 'left':
            self.label_at['hinge'] = QtCore.QPointF(CX - 118, ARM_BASE_Y + 26)
            self.label_at['fc'] = QtCore.QPointF(CX - 118, 300)
        elif self.label_side == 'right':
            self.label_at['hinge'] = QtCore.QPointF(CX + 66, ARM_BASE_Y + 26)
            self.label_at['fc'] = QtCore.QPointF(CX + 66, 300)
        else:
            self.label_at['hinge'] = QtCore.QPointF(70, 400)
            self.label_at['fc'] = QtCore.QPointF(70, 450)

    def _transform(self):
        scale = min(self.width() / SIZE.width(),
                    self.height() / SIZE.height())
        t = QtGui.QTransform()
        t.translate((self.width() - SIZE.width() * scale) / 2,
                    (self.height() - SIZE.height() * scale) / 2)
        t.scale(scale, scale)
        return t

    def part_at(self, pos):
        inverse, ok = self._transform().inverted()
        if not ok:
            return ''
        point = inverse.map(QtCore.QPointF(pos))
        if self.hinge_area.contains(point):
            return 'hinge'
        for part, _, path in self.boxes:
            if path.contains(point):
                return part
        for part, path, _ in self.sites:
            if path.contains(point):
                return part
        return ''

    # state

    def set_fusions(self, kinds):
        '''
        {fusion part: unit kind or None}
        '''
        self.fusions = dict(kinds)
        self._build_shapes()
        self.update()

    def set_formats(self, formats, titles, position='fv'):
        self.formats = dict(formats)
        self.position = position
        self.titles = dict(antibody.PART_LABELS, **titles)
        self._build_shapes()
        self.update()

    def set_part(self, part, assigned, text):
        self.assigned[part] = assigned
        self.texts[part] = text
        self.update()

    # painting

    def _fill(self, part, kind):
        if self.assigned.get(part) or part.startswith('fusion_'):
            key = 'fusion' if part.startswith('fusion_') else part
            color = _color(antibody.PART_COLORS[key],
                           antibody.LIGHT_TINT if kind == 'light' else 0.0)
        else:
            color = QtGui.QColor(theme.COLORS['antibody_part'])
        if part == self.hover:
            color = color.lighter(118)
        return color

    def paintEvent(self, event):
        painter = QtGui.QPainter(self)
        painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing)
        painter.setTransform(self._transform())
        edge = QtGui.QColor(theme.COLORS['antibody_edge'])
        # hinge under the boxes
        hinge = self._fill('hinge', 'heavy') if self.assigned.get(
            'hinge') or self.hover == 'hinge' else \
            QtGui.QColor(theme.COLORS['antibody_part'])
        for width, color in ((5.0, edge), (2.6, hinge)):
            pen = QtGui.QPen(color, width)
            pen.setCapStyle(Qt.PenCapStyle.RoundCap)
            pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
            painter.setPen(pen)
            for path in self.hinge_lines:
                painter.drawPath(path)
        # linkers of fused units
        linker = _color(antibody.PART_COLORS['hinge'])
        for width, color in ((4.0, edge), (2.0, linker)):
            pen = QtGui.QPen(color, width)
            pen.setCapStyle(Qt.PenCapStyle.RoundCap)
            painter.setPen(pen)
            for _, path in self.linker_lines:
                painter.drawPath(path)
        painter.setPen(QtGui.QPen(edge, 2.0))
        for part, kind, path in self.boxes:
            painter.setBrush(self._fill(part, kind))
            painter.drawPath(path)
        # sites to add a fused unit: dashed, "+"
        for part, path, center in self.sites:
            color = QtGui.QColor(theme.COLORS['text_bright' if part ==
                                              self.hover else 'text_dim'])
            pen = QtGui.QPen(color, 1.4, Qt.PenStyle.DashLine)
            painter.setPen(pen)
            painter.setBrush(Qt.BrushStyle.NoBrush)
            painter.drawPath(path)
            painter.drawLine(center - QtCore.QPointF(6, 0),
                             center + QtCore.QPointF(6, 0))
            painter.drawLine(center - QtCore.QPointF(0, 6),
                             center + QtCore.QPointF(0, 6))
        # labels
        font = QtGui.QFont(self.font())
        font.setPointSizeF(theme.FONT_SIZE + 1)
        small = QtGui.QFont(self.font())
        small.setPointSizeF(theme.FONT_SIZE - 1)
        for part, at in self.label_at.items():
            title = self.titles[part]
            text = self.texts.get(part, '')
            color = QtGui.QColor(theme.COLORS['text_bright'] if part ==
                                 self.hover else theme.COLORS['text'])
            box = QtCore.QRectF(at.x() - 105, at.y() - 22, 210, 22)
            note = QtCore.QRectF(at.x() - 105, at.y(), 210, 34)
            align = Qt.AlignmentFlag.AlignHCenter
            if part in ('hinge', 'fc') and self.label_side == 'left':
                # (left of the hinge or Fc, inside the picture: wraps)
                box = QtCore.QRectF(8, at.y() - 22, at.x() + 52 - 8, 22)
                note = QtCore.QRectF(8, at.y(), at.x() + 52 - 8, 48)
                align = Qt.AlignmentFlag.AlignRight
            elif part in ('hinge', 'fc') and self.label_side == 'right':
                width = SIZE.width() - 8 - at.x()
                box = QtCore.QRectF(at.x(), at.y() - 22, width, 22)
                note = QtCore.QRectF(at.x(), at.y(), width, 48)
                align = Qt.AlignmentFlag.AlignLeft
            elif part in ('hinge', 'fc'):
                box = QtCore.QRectF(8, at.y() - 22, 120, 22)
                note = QtCore.QRectF(8, at.y(), 120, 30)
                align = Qt.AlignmentFlag.AlignLeft
            painter.setFont(font)
            painter.setPen(color)
            painter.drawText(box, int(align) |
                             int(Qt.AlignmentFlag.AlignBottom), title)
            painter.setFont(small)
            painter.setPen(QtGui.QColor(theme.COLORS['text_dim']))
            painter.drawText(note, int(align) |
                             int(Qt.AlignmentFlag.AlignTop) |
                             int(Qt.TextFlag.TextWordWrap), text)

    # mouse

    def mouseMoveEvent(self, event):
        pos = event.position() if hasattr(event, 'position') else event.pos()
        part = self.part_at(pos.toPoint())
        if part != self.hover:
            self.hover = part
            self.setCursor(Qt.CursorShape.PointingHandCursor if part else
                           Qt.CursorShape.ArrowCursor)
            if part.startswith('fusion_') and not self.fusions.get(part):
                tip = 'Add an scFv or VHH to the C-terminus of %s' % \
                    antibody.FUSION_SITES[part[len('fusion_'):]]
            elif part:
                tip = '%s: %s' % (self.titles.get(part, TITLES.get(part)),
                                  self.texts.get(part, ''))
            else:
                tip = ''
            self.setToolTip(tip)
            self.update()

    def leaveEvent(self, event):
        self.hover = ''
        self.update()
        super().leaveEvent(event)

    def mousePressEvent(self, event):
        if event.button() != Qt.MouseButton.LeftButton:
            return
        pos = event.position() if hasattr(event, 'position') else event.pos()
        part = self.part_at(pos.toPoint())
        if part:
            self.clicked.emit(part)


class PartDialog(QtWidgets.QDialog):
    '''
    Asks for a part's structure: a loaded object, a file or a PDB ID, and
    its chains (recognized from the sequences)
    '''

    def __init__(self, part, spec, parent=None, _self=cmd, fmt='fab',
                 title='', same_as=None):
        super().__init__(parent)
        self.part = part
        self.fmt = fmt
        self.cmd = _self
        self.title = title or TITLES[part]
        self.setObjectName('antibody_part_dialog')
        self.setWindowTitle('%s structure' % self.title)
        self.setMinimumWidth(480)
        layout = QtWidgets.QVBoxLayout(self)
        if part in FUSIONS:
            fmt = self.fmt = 'fusion'
        note = QtWidgets.QLabel(ARM_DESCRIPTIONS[fmt] if part in ARMS or
                                part in FUSIONS else DESCRIPTIONS[part], self)
        note.setObjectName('prep_note')
        note.setWordWrap(True)
        layout.addWidget(note)

        self.same = None
        if part in ARMS and same_as:
            self.same = QtWidgets.QCheckBox('Same as %s' % same_as, self)
            self.same.setChecked(spec is None)
            self.same.toggled.connect(self._update_enabled)
            layout.addWidget(self.same)
        self.linker = self.both = self.cdrs = None
        self.shape = None
        if part in FUSIONS or (part in ARMS and fmt in ('scfv', 'vhh')):
            # any sequence typed in; presets fill it in
            grid = QtWidgets.QGridLayout()
            grid.addWidget(QtWidgets.QLabel(
                'Linker' if part in FUSIONS else 'Linker to the hinge',
                self), 0, 0)
            self.linker = QtWidgets.QLineEdit(self)
            self.linker.setObjectName('antibody_linker')
            self.linker.setPlaceholderText(
                'type any sequence, e.g. GGGGSGGGGS, (G4S)3, (EAAAK)3; '
                'empty: none')
            self.linker.setToolTip('Residues between the chain\'s end and '
                                   'the unit, built by the refinement'
                                   if part in FUSIONS else
                                   'Residues between the %s and the hinge, '
                                   'built by the refinement (on the '
                                   'template\'s Fv, the hinge alone may be '
                                   'nearly stretched)' % self.title)
            default = antibody.DEFAULT_LINKER if part in FUSIONS else ''
            self.linker.setText(default if spec is None else
                                spec.get('linker', ''))
            grid.addWidget(self.linker, 0, 1, 1, 2)
            self.linker_preset = QtWidgets.QComboBox(self)
            self.linker_preset.addItem('Presets', None)
            for label, text in antibody.LINKERS:
                self.linker_preset.addItem(label, text)
            self.linker_preset.activated.connect(self._use_linker)
            grid.addWidget(self.linker_preset, 1, 1)
            self.linker_note = QtWidgets.QLabel(self)
            self.linker_note.setObjectName('prep_note')
            grid.addWidget(self.linker_note, 1, 2)
            # scFv domain order (an scFv chain: enabled in _object_changed)
            grid.addWidget(QtWidgets.QLabel('scFv domains', self), 2, 0)
            self.orientation = QtWidgets.QComboBox(self)
            self.orientation.setObjectName('antibody_orientation')
            for label, value in ORIENTATION_ITEMS:
                self.orientation.addItem(label, value)
            self.orientation.setToolTip(
                'The order of VH and VL along the chain. Reordering keeps '
                'the same Fv (the domains stay paired as in the structure) '
                'and joins them the other way round, by this linker:')
            self.orientation.setCurrentIndex(self.orientation.findData(
                (spec or {}).get('orientation', '')))
            self.orientation.currentIndexChanged.connect(
                self._orientation_changed)
            grid.addWidget(self.orientation, 2, 1)
            self.scfv_linker = QtWidgets.QLineEdit(self)
            self.scfv_linker.setObjectName('antibody_scfv_linker')
            self.scfv_linker.setPlaceholderText('linker between them')
            self.scfv_linker.setText((spec or {}).get(
                'scfv_linker', antibody.DEFAULT_SCFV_LINKER))
            self.scfv_linker.setToolTip('Residues joining VH and VL in the '
                                        'new order, built by the refinement '
                                        '(the structure\'s own linker can\'t '
                                        'span the other way round)')
            grid.addWidget(self.scfv_linker, 2, 2)
            grid.setColumnStretch(2, 1)
            layout.addLayout(grid)
            self.linker.textChanged.connect(self._linker_changed)
            self._linker_changed()
            self._orientation_changed()
            if part in FUSIONS:
                kind = part[len('fusion_'):-1]
                self.both = QtWidgets.QCheckBox('Both %s chains' % kind,
                                                self)
                self.both.setChecked(spec is None or bool(spec.get('both')))
                layout.addWidget(self.both)
                row = QtWidgets.QHBoxLayout()
                row.addWidget(QtWidgets.QLabel('Linker shape', self))
                self.shape = QtWidgets.QComboBox(self)
                self.shape.setObjectName('antibody_linker_shape')
                for label, value in SHAPE_ITEMS:
                    self.shape.addItem(label, value)
                self.shape.setCurrentIndex(self.shape.findData(
                    (spec or {}).get('linker_shape', 'relaxed')))
                self.shape.setToolTip(
                    'Where the unit sits: near the chain\'s end with its '
                    'linker coiled, at the length of the stretched linker, '
                    'or stretched out and then drawn in by a short '
                    'dynamics run, which lets the linker coil itself '
                    '(slower)')
                row.addWidget(self.shape, 1)
                layout.addLayout(row)
                self.cdrs = QtWidgets.QCheckBox(
                    'Binding site facing outwards', self)
                self.cdrs.setObjectName('antibody_cdrs')
                self.cdrs.setChecked((spec or {}).get('cdrs', 'out') ==
                                     'out')
                self.cdrs.setToolTip(
                    'Turn the unit so its CDRs face the solvent rather than '
                    'the antibody, as far as that doesn\'t make it clash')
                layout.addWidget(self.cdrs)
        self.modes = None
        if part == 'hinge':
            self.modes = QtWidgets.QButtonGroup(self)
            row = QtWidgets.QHBoxLayout()
            for k, label in enumerate(('From the template', 'Sequence',
                                       'Structure')):
                button = QtWidgets.QRadioButton(label, self)
                self.modes.addButton(button, k)
                row.addWidget(button)
            row.addStretch(1)
            layout.addLayout(row)
            self.sequence_box = QtWidgets.QWidget(self)
            seq_row = QtWidgets.QHBoxLayout(self.sequence_box)
            seq_row.setContentsMargins(0, 0, 0, 0)
            self.preset = QtWidgets.QComboBox(self.sequence_box)
            self.preset.addItem('Isotype...', '')
            for label, text in antibody.HINGES:
                self.preset.addItem(label, text)
            self.preset.activated.connect(self._use_preset)
            seq_row.addWidget(self.preset)
            self.sequence = QtWidgets.QLineEdit(self.sequence_box)
            self.sequence.setPlaceholderText('e.g. EPKSCDKTHTCPPCP')
            self.sequence.setToolTip('One-letter codes; residues the Fab or '
                                     'the Fc already have are left out')
            seq_row.addWidget(self.sequence, 1)
            layout.addWidget(self.sequence_box)

        # structure
        self.structure_box = QtWidgets.QGroupBox('Structure', self)
        self.structure_box.setObjectName('prep_step')
        grid = QtWidgets.QGridLayout(self.structure_box)
        self.object = QtWidgets.QComboBox(self.structure_box)
        self.object.setMinimumContentsLength(14)
        self.object.currentIndexChanged.connect(self._object_changed)
        grid.addWidget(QtWidgets.QLabel('Object', self.structure_box), 0, 0)
        grid.addWidget(self.object, 0, 1)
        open_button = QtWidgets.QPushButton('Open File...',
                                            self.structure_box)
        open_button.clicked.connect(self._open_file)
        grid.addWidget(open_button, 0, 2)
        grid.addWidget(QtWidgets.QLabel('PDB ID', self.structure_box), 1, 0)
        self.pdb_id = QtWidgets.QLineEdit(self.structure_box)
        self.pdb_id.setPlaceholderText('e.g. 1hzh')
        self.pdb_id.setMaxLength(4)
        self.pdb_id.returnPressed.connect(self._fetch)
        grid.addWidget(self.pdb_id, 1, 1)
        fetch = QtWidgets.QPushButton('Fetch', self.structure_box)
        fetch.clicked.connect(self._fetch)
        grid.addWidget(fetch, 1, 2)
        if part == 'template':
            self.known = QtWidgets.QComboBox(self.structure_box)
            self.known.addItem('Standard templates...', '')
            for code, label in antibody.TEMPLATES:
                self.known.addItem(label, code)
            self.known.activated.connect(self._fetch_known)
            grid.addWidget(QtWidgets.QLabel('Or', self.structure_box), 2, 0)
            grid.addWidget(self.known, 2, 1, 1, 2)
        self.chains = []
        row = grid.rowCount()
        for k, label in enumerate(ARM_CHAIN_LABELS[fmt]
                                  if part in ARMS or part in FUSIONS
                                  else CHAIN_LABELS[part]):
            combo = QtWidgets.QComboBox(self.structure_box)
            self.chains.append(combo)
            grid.addWidget(QtWidgets.QLabel(label, self.structure_box),
                           row + k, 0)
            grid.addWidget(combo, row + k, 1, 1, 2)
        grid.setColumnStretch(1, 1)
        self.chain_note = QtWidgets.QLabel(self.structure_box)
        self.chain_note.setObjectName('prep_note')
        self.chain_note.setWordWrap(True)
        grid.addWidget(self.chain_note, grid.rowCount(), 0, 1, 3)
        layout.addWidget(self.structure_box)

        buttons = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Ok |
            QtWidgets.QDialogButtonBox.StandardButton.Cancel, self)
        buttons.accepted.connect(self._accept)
        buttons.rejected.connect(self.reject)
        if part != 'template' and self.same is None:
            clear = buttons.addButton(
                'Clear', QtWidgets.QDialogButtonBox.ButtonRole.ResetRole)
            clear.clicked.connect(self._clear)
        self.cleared = False
        layout.addWidget(buttons)

        self._fill_objects(spec['object'] if spec and spec.get('object')
                           else '')
        if spec and spec.get('chains'):
            for combo, chain in zip(self.chains, spec['chains']):
                index = combo.findData(chain)
                if index >= 0:
                    combo.setCurrentIndex(index)
        if part == 'template' and spec and spec.get('object') and \
                self.object.findText(spec['object']) < 0:
            self.pdb_id.setText(spec['object'])
        if self.modes is not None:
            mode = (spec or {}).get('mode', 'template')
            self.modes.button(('template', 'sequence', 'structure').index(
                mode)).setChecked(True)
            self.sequence.setText((spec or {}).get('sequence', ''))
            self.modes.idToggled.connect(self._update_enabled)
        self._update_enabled()

    # structure

    def _fill_objects(self, current=''):
        names = [n for n in self.cmd.get_names('objects') or []
                 if self.cmd.get_type(n) == 'object:molecule' and
                 not n.startswith('_')]
        self.object.blockSignals(True)
        self.object.clear()
        self.object.addItems(names)
        if current in names:
            self.object.setCurrentText(current)
        self.object.blockSignals(False)
        self._object_changed()

    def _object_changed(self, *_):
        obj = self.object.currentText()
        for combo in self.chains:
            combo.clear()
        self.chain_note.setText('')
        if not obj:
            return
        chains = antibody.protein_chains(obj, _self=self.cmd)
        kinds = {'heavy': 'heavy', 'light': 'light', 'fc': 'Fc',
                 'scfv': 'scFv', 'vhh': 'VHH'}
        for combo in self.chains:
            for chain, (sequence, kind) in chains.items():
                combo.addItem('%s  (%s%d aa)' % (
                    chain or '-', kinds[kind] + ', ' if kind else '',
                    len(sequence)), chain)
        try:
            if self.part in FUSIONS:
                found = (antibody.detect_single(obj, 'scfv', self.cmd) or
                         antibody.detect_single(obj, 'vhh', self.cmd),)
            elif self.part in ARMS and self.fmt != 'fab':
                found = (antibody.detect_single(obj, self.fmt, self.cmd),)
            else:
                found = {'fab1': antibody.detect_fab,
                         'fab2': antibody.detect_fab,
                         'fc': antibody.detect_fc,
                         'hinge': antibody.detect_hinge,
                         'template': antibody.detect_template}[self.part](
                             obj, _self=self.cmd)
        except ValueError as e:
            self.chain_note.setText(str(e))
            return
        for combo, chain in zip(self.chains, found):
            index = combo.findData(chain)
            if index >= 0:
                combo.setCurrentIndex(index)
        if not all(found):
            self.chain_note.setText('Some chains weren\'t recognized: '
                                    'choose them')
        self._update_orientation(chains)

    def _update_orientation(self, chains=None):
        '''
        The domain order is an scFv's; a VHH has one domain
        '''
        if getattr(self, 'orientation', None) is None:
            return
        chain = self.chains[0].currentData() if self.chains else None
        if chains is None:
            chains = antibody.protein_chains(self.object.currentText(),
                                             _self=self.cmd) \
                if self.object.currentText() else {}
        scfv = chains.get(chain, ('', ''))[1] == 'scfv'
        self.orientation.setEnabled(scfv)
        self._orientation_changed()
        if not scfv:
            self.orientation.setCurrentIndex(
                self.orientation.findData(''))

    def _load(self, action):
        before = set(self.cmd.get_names('objects') or [])
        QtWidgets.QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        try:
            action()
        except Exception as e:
            QtWidgets.QApplication.restoreOverrideCursor()
            self.chain_note.setText('Failed: %s' % e)
            return
        QtWidgets.QApplication.restoreOverrideCursor()
        new = [n for n in self.cmd.get_names('objects') or []
               if n not in before]
        self._fill_objects(new[0] if new else self.object.currentText())

    def _open_file(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, 'Open %s structure' % self.title, '',
            'Structures (*.pdb *.cif *.mmcif *.ent *.pdb.gz *.cif.gz '
            '*.bcif);;All files (*)')
        if path:
            self._load(lambda: self.cmd.load(path, quiet=1))

    def _fetch(self):
        code = self.pdb_id.text().strip().lower()
        if len(code) == 4:
            if code in (self.cmd.get_names('objects') or []):
                self._fill_objects(code)
            else:
                self._load(lambda: self.cmd.fetch(code, quiet=1))

    def _fetch_known(self, index):
        code = self.known.itemData(index)
        if code:
            self.pdb_id.setText(code)
            self._fetch()
        self.known.setCurrentIndex(0)

    def _linker_changed(self, *_):
        '''
        The linker typed in: its residues, or what's wrong
        '''
        try:
            text = antibody.parse_linker(self.linker.text())
        except ValueError as e:
            self.linker_note.setText(str(e).split(': ', 1)[-1])
            return
        if not text:
            self.linker_note.setText('no linker (joined directly)')
        else:
            self.linker_note.setText('%d residues%s' % (
                len(text), ': ' + text if text != self.linker.text().upper()
                .replace(' ', '') else ''))

    def _orientation_changed(self, *_):
        '''
        The linker between the domains matters only when reordering
        '''
        self.scfv_linker.setEnabled(bool(self.orientation.currentData()))

    def _use_linker(self, index):
        text = self.linker_preset.itemData(index)
        if text is not None:
            self.linker.setText(text)
        self.linker_preset.setCurrentIndex(0)

    def _use_preset(self, index):
        text = self.preset.itemData(index)
        if text:
            self.sequence.setText(text)
            self.modes.button(1).setChecked(True)
        self.preset.setCurrentIndex(0)

    def _mode(self):
        if self.modes is None:
            return ''
        return ('template', 'sequence', 'structure')[self.modes.checkedId()]

    def _update_enabled(self, *_):
        on = True
        if self.same is not None:
            on = not self.same.isChecked()
        if self.modes is not None:
            on = self._mode() == 'structure'
            self.sequence_box.setEnabled(self._mode() == 'sequence')
        self.structure_box.setEnabled(on)

    # result

    def _clear(self):
        self.cleared = True
        self.accept()

    def _accept(self):
        if self.structure_box.isEnabled():
            if not self.object.currentText():
                self.chain_note.setText('Choose an object, open a file or '
                                        'fetch a PDB entry')
                return
            if not all(c.currentData() is not None for c in self.chains):
                self.chain_note.setText('Choose the chains')
                return
        if self._mode() == 'sequence' and not self.sequence.text().strip():
            self.sequence.setFocus()
            return
        if self.linker is not None:
            try:
                antibody.parse_linker(self.linker.text())
                antibody.parse_linker(self.scfv_linker.text())
            except ValueError as e:
                self.chain_note.setText(str(e))
                return
        self.accept()

    def spec(self):
        '''
        None (cleared, arm 2 same as arm 1) or {'object', 'chains'} (and for
        the hinge 'mode', 'sequence')
        '''
        if self.cleared or (self.same is not None and self.same.isChecked()):
            return None
        result = {'object': self.object.currentText(),
                  'chains': tuple(c.currentData() for c in self.chains)}
        if self.modes is not None:
            result['mode'] = self._mode()
            result['sequence'] = ''.join(
                ch for ch in self.sequence.text().upper() if ch.isalpha())
            if result['mode'] != 'structure':
                result['object'], result['chains'] = '', ()
        if self.linker is not None:
            result['linker'] = antibody.parse_linker(self.linker.text())
            result['orientation'] = self.orientation.currentData() \
                if self.orientation.isEnabled() else ''
            result['scfv_linker'] = antibody.parse_linker(
                self.scfv_linker.text())
        if self.cdrs is not None:
            result['cdrs'] = 'out' if self.cdrs.isChecked() else 'any'
            result['linker_shape'] = self.shape.currentData()
        if self.both is not None:
            result['both'] = self.both.isChecked()
            chains = antibody.protein_chains(result['object'], self.cmd)
            kind = chains.get(result['chains'][0], ('', ''))[1]
            result['kind'] = 'scfv' if kind == 'scfv' else 'vhh'
        return result


def _part_name(assembly, chain, part):
    if part == 'hinge':
        return 'hinge (%s)' % chain
    if part == 'fc':
        return 'Fc (%s)' % chain
    if assembly.formats.get(part) != 'fab':
        return '%s (%s)' % (assembly.labels[part], chain)
    kind = 'light' if chain in assembly.light else 'heavy'
    return '%s %s' % (assembly.labels[part], kind)


class DisulfideDialog(QtWidgets.QDialog):
    '''
    The model's cysteine pairs to choose disulfide bonds from (checked:
    the current choice, or those close enough to form by themselves)
    '''

    def __init__(self, assembly, chosen=None, parent=None):
        super().__init__(parent)
        self.setObjectName('antibody_disulfide_dialog')
        self.setWindowTitle('Disulfide Bonds')
        self.resize(560, 520)
        layout = QtWidgets.QVBoxLayout(self)
        note = QtWidgets.QLabel(
            'Checked pairs become disulfide bonds, the others free '
            'cysteines. Cysteines apart are brought together by the '
            'refinement (far apart, it can strain the hinge). Numbers are '
            'the model\'s.', self)
        note.setObjectName('prep_note')
        note.setWordWrap(True)
        layout.addWidget(note)
        self.pairs = antibody.cysteine_pairs(assembly)
        self.table = QtWidgets.QTableWidget(len(self.pairs), 3, self)
        self.table.setObjectName('antibody_disulfides')
        self.table.setHorizontalHeaderLabels(['Cysteines', 'Parts',
                                              'SG-SG (Å)'])
        self.table.verticalHeader().setVisible(False)
        self.table.setEditTriggers(
            QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionMode(
            QtWidgets.QAbstractItemView.SelectionMode.NoSelection)
        self.table.horizontalHeader().setStretchLastSection(True)
        known = set()
        chosen_set = None if chosen is None else \
            {frozenset(pair) for pair in chosen}
        for row, pair in enumerate(self.pairs):
            key = frozenset((pair['a'], pair['b']))
            known.add(key)
            item = QtWidgets.QTableWidgetItem(
                antibody.pair_label((pair['a'], pair['b'])))
            item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            on = pair['bonded'] if chosen_set is None else key in chosen_set
            item.setCheckState(Qt.CheckState.Checked if on else
                               Qt.CheckState.Unchecked)
            self.table.setItem(row, 0, item)
            self.table.setItem(row, 1, QtWidgets.QTableWidgetItem(
                ' - '.join(_part_name(assembly, c[0], p) for c, p in
                           zip((pair['a'], pair['b']), pair['parts']))))
            d = pair['distance']
            item = QtWidgets.QTableWidgetItem(
                'to build' if d is None else '%.1f' % d)
            item.setTextAlignment(int(Qt.AlignmentFlag.AlignRight) |
                                  int(Qt.AlignmentFlag.AlignVCenter))
            self.table.setItem(row, 2, item)
        self.table.resizeColumnsToContents()
        layout.addWidget(self.table, 1)
        row = QtWidgets.QHBoxLayout()
        row.addWidget(QtWidgets.QLabel('Other pairs', self))
        self.other = QtWidgets.QLineEdit(self)
        self.other.setPlaceholderText('e.g. H/232-K/232')
        if chosen:
            self.other.setText(', '.join(
                antibody.pair_label(pair) for pair in chosen
                if frozenset(pair) not in known))
        row.addWidget(self.other, 1)
        layout.addLayout(row)
        self.error = QtWidgets.QLabel(self)
        self.error.setObjectName('prep_note')
        layout.addWidget(self.error)
        buttons = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Ok |
            QtWidgets.QDialogButtonBox.StandardButton.Cancel, self)
        automatic = buttons.addButton(
            'Automatic', QtWidgets.QDialogButtonBox.ButtonRole.ResetRole)
        automatic.setToolTip('Those close enough to form by themselves')
        automatic.clicked.connect(self._automatic)
        buttons.accepted.connect(self._accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)
        self.automatic = False
        self.chosen = chosen

    def _automatic(self):
        self.automatic = True
        self.chosen = None
        self.accept()

    def _accept(self):
        from pymolx.prep import parse_disulfides
        try:
            other = parse_disulfides(self.other.text()) or []
        except ValueError as e:
            self.error.setText(str(e))
            return
        chosen = [(pair['a'], pair['b']) for row, pair in
                  enumerate(self.pairs)
                  if self.table.item(row, 0).checkState() ==
                  Qt.CheckState.Checked]
        self.chosen = chosen + [p for p in other if p not in chosen]
        self.accept()


class RefineThread(QtCore.QThread):
    '''
    Runs the refinement (PrepJob.run, no PyMOL calls) off the GUI thread
    '''

    def __init__(self, job, parent=None):
        super().__init__(parent)
        self.job = job
        self.error = None
        self.cancelled = False

    def run(self):
        try:
            self.job.run()
        except prep.Cancelled:
            self.cancelled = True
        except Exception as e:
            self.error = e
            traceback.print_exc()


class AntibodyDialog(QtWidgets.QDialog):

    def __init__(self, parent=None, _self=cmd):
        super().__init__(parent)
        self.cmd = _self
        self.specs = {'fab1': None, 'fab2': None, 'hinge': None, 'fc': None,
                      'template': {'object': antibody.DEFAULT_TEMPLATE,
                                   'chains': ()}}
        self.specs.update({part: None for part in FUSIONS})
        self.formats = {'fab1': 'fab', 'fab2': 'fab'}
        self.thread = None
        self.assembly = None
        self.setObjectName('antibody_dialog')
        self.setWindowTitle('Antibody Builder')
        self.resize(940, 640)

        layout = QtWidgets.QHBoxLayout(self)
        left = QtWidgets.QVBoxLayout()
        layout.addLayout(left, 3)
        note = QtWidgets.QLabel('Click a part of the antibody to give its '
                                'structure', self)
        note.setObjectName('prep_note')
        left.addWidget(note)
        self.schematic = Schematic(self)
        self.schematic.clicked.connect(self.edit_part)
        left.addWidget(self.schematic, 1)

        right = QtWidgets.QVBoxLayout()
        layout.addLayout(right, 2)
        box = QtWidgets.QGroupBox('Format', self)
        box.setObjectName('prep_step')
        grid = QtWidgets.QGridLayout(box)
        row = QtWidgets.QHBoxLayout()
        grid.addLayout(row, 0, 0, 1, 2)
        self.format_boxes = {}
        for k, part in enumerate(ARMS):
            row.addWidget(QtWidgets.QLabel('Arm %d' % (k + 1), box))
            combo = QtWidgets.QComboBox(box)
            combo.setObjectName('antibody_%s_format' % part)
            for label, value in FORMAT_ITEMS:
                combo.addItem(label, value)
            combo.setToolTip('Arm %d (%s in the picture): a Fab, an scFv, '
                             'a single domain (VHH) or none (one-armed)' % (
                                 k + 1, 'left' if k == 0 else 'right'))
            combo.currentIndexChanged.connect(
                lambda _=0, p=part: self._format_changed(p))
            self.format_boxes[part] = combo
            row.addWidget(combo, 1)
        grid.addWidget(QtWidgets.QLabel('scFv, VHH', box), 1, 0)
        self.position = QtWidgets.QComboBox(box)
        self.position.setObjectName('antibody_single_position')
        for label, value in POSITIONS:
            self.position.addItem(label, value)
        self.position.setToolTip(
            'Where scFv and VHH arms go: superposed on the template arm\'s '
            'Fv, or then moved next to the hinge into the place of CH1/CL '
            '(as in scFv-Fc)')
        self.position.currentIndexChanged.connect(self._show_specs)
        grid.addWidget(self.position, 1, 1)
        grid.setColumnStretch(1, 1)
        right.addWidget(box)

        box = QtWidgets.QGroupBox('Parts', self)
        box.setObjectName('prep_step')
        grid = QtWidgets.QGridLayout(box)
        self.part_labels = {}
        self.part_titles = {}
        self.part_buttons = {}
        for row, part in enumerate(('fab1', 'fab2', 'hinge', 'fc',
                                    'template', 'fusions')):
            title = QtWidgets.QLabel(TITLES[part], box)
            self.part_titles[part] = title
            grid.addWidget(title, row, 0)
            label = QtWidgets.QLabel(box)
            label.setObjectName('prep_note')
            self.part_labels[part] = label
            grid.addWidget(label, row, 1)
            button = QtWidgets.QPushButton('Change...', box)
            button.setObjectName('antibody_%s_button' % part)
            if part == 'fusions':
                button.setText('Add...')
                button.setToolTip('An scFv or VHH appended to a chain\'s '
                                  'C-terminus (or click the dashed boxes '
                                  'under the Fc)')
                self.fusion_menu = QtWidgets.QMenu(button)
                button.setMenu(self.fusion_menu)
                self.fusion_menu.aboutToShow.connect(self._fill_fusion_menu)
            else:
                button.clicked.connect(
                    lambda _=False, p=part: self.edit_part(p))
            self.part_buttons[part] = button
            grid.addWidget(button, row, 2)
        grid.setColumnStretch(1, 1)
        right.addWidget(box)

        box = QtWidgets.QGroupBox('Model', self)
        box.setObjectName('prep_step')
        grid = QtWidgets.QGridLayout(box)
        self.name = QtWidgets.QLineEdit('igg_model', box)
        grid.addWidget(QtWidgets.QLabel('Name', box), 0, 0)
        grid.addWidget(self.name, 0, 1)
        self.chain_ids = QtWidgets.QLineEdit('HLKM', box)
        self.chain_ids.setToolTip('Chain IDs of the model: arm 1 heavy (or '
                                  'scFv, VHH), arm 1 light, arm 2 heavy, '
                                  'arm 2 light')
        grid.addWidget(QtWidgets.QLabel('Chains', box), 1, 0)
        grid.addWidget(self.chain_ids, 1, 1)
        self.glycans = QtWidgets.QCheckBox('Keep glycans', box)
        self.glycans.setChecked(True)
        self.glycans.setToolTip('Glycans and other groups bonded to the '
                                'chains go with them')
        grid.addWidget(self.glycans, 2, 0, 1, 2)
        self.refine = QtWidgets.QCheckBox('Build missing residues and '
                                          'minimize (OpenMM)', box)
        self.refine.setChecked(True)
        self.refine.setToolTip('Missing residues and side chains are built, '
                               'hydrogens added and the model minimized; '
                               'the hinge and the junctions move freely, '
                               'the rest is restrained')
        grid.addWidget(self.refine, 3, 0, 1, 2)
        grid.addWidget(QtWidgets.QLabel('Disulfides', box), 4, 0)
        row = QtWidgets.QHBoxLayout()
        self.disulfide_label = QtWidgets.QLabel(box)
        self.disulfide_label.setObjectName('prep_note')
        row.addWidget(self.disulfide_label, 1)
        self.disulfide_button = QtWidgets.QPushButton('Choose...', box)
        self.disulfide_button.setToolTip('Which cysteine pairs form '
                                         'disulfide bonds')
        self.disulfide_button.clicked.connect(self.choose_disulfides)
        row.addWidget(self.disulfide_button)
        grid.addLayout(row, 4, 1)
        right.addWidget(box)
        self.disulfides = None        # chosen pairs, or None: automatic
        self.disulfide_args = None    # for the parts they were chosen for

        row = QtWidgets.QHBoxLayout()
        self.status = QtWidgets.QLabel(self)
        self.status.setObjectName('prep_note')
        self.status.setWordWrap(True)
        row.addWidget(self.status, 1)
        self.cancel_button = QtWidgets.QPushButton('Cancel', self)
        self.cancel_button.setEnabled(False)
        self.cancel_button.clicked.connect(self.cancel)
        row.addWidget(self.cancel_button)
        self.build_button = QtWidgets.QPushButton('Build Antibody', self)
        self.build_button.setDefault(True)
        self.build_button.clicked.connect(self.build)
        row.addWidget(self.build_button)
        right.addLayout(row)
        self.log = QtWidgets.QPlainTextEdit(self)
        self.log.setObjectName('antibody_log')
        self.log.setReadOnly(True)
        self.log.setFont(theme.console_font())
        right.addWidget(self.log, 1)

        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self._poll)
        self._log_lines = 0
        self._show_specs()

    # parts

    def _show_disulfides(self):
        self.disulfide_label.setText(
            'by distance (automatic)' if self.disulfides is None else
            '%d chosen' % len(self.disulfides))

    def choose_disulfides(self):
        if self.thread is not None:
            return
        QtWidgets.QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        try:
            args = self._arguments()
            assembly = antibody.assemble(_self=self.cmd, **args)
        except Exception as e:
            self.status.setText('Disulfides: %s' % e)
            return
        finally:
            QtWidgets.QApplication.restoreOverrideCursor()
        self._check_objects()   # (the template may have been fetched)
        chosen = self.disulfides if self.disulfide_args == args else None
        dialog = DisulfideDialog(assembly, chosen, self)
        if dialog.exec() == QtWidgets.QDialog.DialogCode.Accepted:
            self.disulfides = dialog.chosen
            self.disulfide_args = None if dialog.chosen is None else args
            self._show_disulfides()

    def _labels(self):
        return antibody.arm_labels([self.formats[p] for p in ARMS])

    def _same_as(self, part):
        '''
        For arm 2 of the same format as arm 1: arm 1's label
        '''
        if part == 'fab2' and self.formats['fab1'] == self.formats['fab2']:
            return self._labels()['fab1']
        return None

    def _format_changed(self, part):
        fmt = self.format_boxes[part].currentData()
        other = self.formats['fab2' if part == 'fab1' else 'fab1']
        if fmt == 'none' and other == 'none':
            self.status.setText('At least one arm')
            self.format_boxes[part].setCurrentIndex(
                self.format_boxes[part].findData(self.formats[part]))
            return
        if fmt != self.formats[part]:
            self.specs[part] = None     # (another kind of structure)
            self.formats[part] = fmt
        self._show_specs()

    def _fill_fusion_menu(self):
        self.fusion_menu.clear()
        for part in FUSIONS:
            site = part[len('fusion_'):]
            spec = self.specs[part]
            action = self.fusion_menu.addAction(
                '%s%s' % (antibody.FUSION_SITES[site].capitalize(),
                          '  (%s)' % summary(part, spec) if spec else ''))
            k = int(site[-1])
            action.setEnabled(site.startswith('heavy') or
                              self.formats['fab%d' % k] == 'fab')
            action.triggered.connect(
                lambda _=False, p=part: self.edit_part(p))

    def _show_fusions(self):
        given = [p for p in FUSIONS if self.specs[p]]
        self.part_labels['fusions'].setText(', '.join(
            '%s on %s' % (antibody.FORMAT_LABELS[self.specs[p]['kind']],
                          antibody.FUSION_SITES[p[len('fusion_'):]])
            for p in given) or 'none')
        self.schematic.set_fusions({p: self.specs[p]['kind'] for p in given})
        for part in given:
            self.schematic.set_part(part, True, summary(part,
                                                        self.specs[part]))

    def _show_specs(self):
        if hasattr(self, 'disulfide_label'):
            self._show_disulfides()
        labels = self._labels()
        self.schematic.set_formats(self.formats, labels,
                                   self.position.currentData())
        self.position.setEnabled(any(self.formats[p] in ('scfv', 'vhh')
                                     for p in ARMS))
        # (light chain fusions need a Fab on their arm)
        for k, part in enumerate(ARMS):
            if self.formats[part] != 'fab':
                self.specs['fusion_light%d' % (k + 1)] = None
        self._show_fusions()
        for part, spec in self.specs.items():
            if part in FUSIONS:
                continue
            fmt = self.formats.get(part, '')
            same_as = self._same_as(part)
            text = summary(part, spec, fmt, same_as)
            self.part_labels[part].setText(text)
            if part in ARMS:
                self.part_titles[part].setText(labels[part])
                self.part_buttons[part].setEnabled(fmt != 'none')
            if part != 'template':
                assigned = spec is not None and (
                    part != 'hinge' or spec['mode'] != 'template')
                if part in ARMS and spec is None and same_as:
                    assigned = self.specs['fab1'] is not None
                self.schematic.set_part(part, assigned, text)

    def _check_objects(self):
        '''
        Parts whose objects were deleted are given again
        '''
        names = set(self.cmd.get_names('objects') or [])
        for part, spec in self.specs.items():
            if spec and spec.get('object') and spec['object'] not in names \
                    and not (part == 'template' and not spec.get('chains')):
                self.specs[part] = None if part != 'template' else {
                    'object': antibody.DEFAULT_TEMPLATE, 'chains': ()}
        self._show_specs()

    def edit_part(self, part):
        if self.thread is not None or self.formats.get(part) == 'none':
            return
        dialog = PartDialog(part, self.specs[part], self, _self=self.cmd,
                            fmt=self.formats.get(part, 'fab'),
                            title=self._labels().get(part, ''),
                            same_as=self._same_as(part))
        if dialog.exec() == QtWidgets.QDialog.DialogCode.Accepted:
            spec = dialog.spec()
            if part == 'template' and spec is None:
                return
            if part == 'hinge' and spec and spec['mode'] == 'template':
                spec = None
            self.specs[part] = spec
            if part in FUSIONS and spec and spec.get('both'):
                other = part[:-1] + ('2' if part.endswith('1') else '1')
                if other in self.specs and (
                        other.startswith('fusion_heavy') or self.formats[
                            'fab' + other[-1]] == 'fab'):
                    self.specs[other] = dict(spec)
            self._show_specs()

    def showEvent(self, event):
        self._check_objects()
        super().showEvent(event)

    # building

    def _arguments(self):
        s = self.specs
        labels = self._labels()
        args = {'arms': '+'.join(self.formats[p] for p in ARMS),
                'single_position': self.position.currentData()}
        linkers = ['', '']
        orientations = ['', '']
        scfv_linkers = [antibody.DEFAULT_SCFV_LINKER] * 2
        for part in ARMS:
            if self.formats[part] == 'none':
                continue
            if s[part] is None:
                if self._same_as(part) and s['fab1'] is not None:
                    # (arm 1's structure again)
                    linkers[1] = s['fab1'].get('linker', '')
                    orientations[1] = s['fab1'].get('orientation', '')
                    scfv_linkers[1] = s['fab1'].get(
                        'scfv_linker', antibody.DEFAULT_SCFV_LINKER)
                    continue
                raise ValueError('give the %s structure (click it)' %
                                 labels[part])
            args[part] = s[part]['object']
            args[part + '_chains'] = ','.join(s[part]['chains'])
            linkers[ARMS.index(part)] = s[part].get('linker', '')
            orientations[ARMS.index(part)] = s[part].get('orientation', '')
            scfv_linkers[ARMS.index(part)] = s[part].get(
                'scfv_linker', antibody.DEFAULT_SCFV_LINKER)
        args['arm_linkers'] = tuple(linkers)
        args['arm_orientations'] = tuple(orientations)
        args['arm_scfv_linkers'] = tuple(scfv_linkers)
        if s['fc'] is None:
            raise ValueError('give the Fc structure (click Fc)')
        args.update({
            'fc': s['fc']['object'],
            'fc_chains': ','.join(s['fc']['chains']),
            'template': s['template']['object'],
            'template_chains': ','.join(s['template']['chains'] or ()),
            'chains': self.chain_ids.text().strip() or 'HLKM',
            'glycans': self.glycans.isChecked()})
        args['fusions'] = [
            {'site': part[len('fusion_'):], 'object': s[part]['object'],
             'chain': s[part]['chains'][0], 'linker': s[part]['linker'],
             'orientation': s[part].get('orientation', ''),
             'scfv_linker': s[part].get('scfv_linker',
                                        antibody.DEFAULT_SCFV_LINKER),
             'cdrs': s[part].get('cdrs', 'out'),
             'linker_shape': s[part].get('linker_shape', 'relaxed')}
            for part in FUSIONS if s[part]]
        hinge = s['hinge']
        if hinge and hinge['mode'] == 'sequence':
            args['hinge_sequence'] = hinge['sequence']
        elif hinge and hinge['mode'] == 'structure':
            args['hinge'] = hinge['object']
            args['hinge_chains'] = ','.join(hinge['chains'])
        return args

    def _say(self, text):
        self.log.appendPlainText(text)

    def build(self):
        self.log.clear()
        self._log_lines = 0
        name = antibody.object_name(self.name.text(), self.cmd)
        QtWidgets.QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        try:
            args = self._arguments()
            self.assembly = antibody.assemble(_self=self.cmd, **args)
        except Exception as e:
            self.status.setText('Failed')
            self._say('Building failed: %s' % e)
            return
        finally:
            QtWidgets.QApplication.restoreOverrideCursor()
        self._check_objects()   # (the template may have been fetched)
        if self.disulfides is not None and self.disulfide_args != args:
            self._say('The parts changed since the disulfide bonds were '
                      'chosen: by distance instead')
            self.disulfides = self.disulfide_args = None
            self._show_disulfides()
        self.assembly.disulfides = self.disulfides
        for line in self.assembly.log:
            self._say(line)
        if self.disulfides is not None:
            self._say('Disulfide bonds chosen: %s' % (', '.join(
                antibody.pair_label(p) for p in self.disulfides) or 'none'))
        if not self.refine.isChecked():
            name = antibody.load(self.assembly, name, self.cmd)
            if self.disulfides is not None:
                antibody.set_disulfides(name, self.disulfides, self.cmd)
            antibody.present(self.assembly, self.cmd)
            self._done()
            return
        self.source = self.cmd.get_unused_name('_ab_assembly')
        antibody.load(self.assembly, self.source, self.cmd)
        try:
            job = antibody.refine_job(self.assembly, self.source, name,
                                      self.cmd)
        except Exception as e:
            self.cmd.delete(self.source)
            self.status.setText('Failed')
            self._say('Refinement failed: %s' % e)
            return
        self._say('Refining (OpenMM)...')
        self.thread = RefineThread(job, self)
        self.thread.finished.connect(self._finished)
        self._set_running(True)
        self.status.setText('Starting...')
        self.thread.start()
        self.timer.start(200)

    def _flush(self, job):
        for line in job.log[self._log_lines:]:
            self._say(line)
        self._log_lines = len(job.log)

    def _poll(self):
        if self.thread is not None:
            self._flush(self.thread.job)
            self.status.setText(self.thread.job.status)

    def _finished(self):
        self.timer.stop()
        thread, self.thread = self.thread, None
        self._set_running(False)
        job = thread.job
        self._flush(job)
        if thread.cancelled or thread.error is not None:
            self.cmd.delete(self.source)
            self.status.setText('Cancelled, nothing changed' if
                                thread.cancelled else 'Failed')
            if thread.error is not None:
                self._say('Refinement failed: %s' % thread.error)
        else:
            count = len(self.assembly.log)
            try:
                antibody.finish(self.assembly, job, self.source, self.cmd)
            except Exception as e:
                self.cmd.delete(self.source)
                self.status.setText('Failed')
                self._say('Loading the model failed: %s' % e)
            else:
                # (finish appends the refinement's log, already shown)
                for line in self.assembly.log[count + len(job.log):]:
                    self._say(line)
                self._done()
        thread.deleteLater()

    def _done(self):
        a = self.assembly
        self.status.setText('Done: "%s"%s' % (
            a.name, ', %d warnings' % len(a.warnings) if a.warnings else ''))
        tags = a.tags()
        present = set(a.parts().values())
        self._say('New model "%s": %s' % (a.name, ', '.join(
            '%s %s' % (a.labels.get(p, p), tag) for p, tag in tags.items()
            if p in present)) + ' (atom property "custom", e.g. select '
                  'custom HNGE)')

    def _set_running(self, running):
        self.build_button.setEnabled(not running)
        self.cancel_button.setEnabled(running)

    def cancel(self):
        if self.thread is not None:
            self.thread.job.cancel()
            self.status.setText('Cancelling...')

    def closeEvent(self, event):
        if self.thread is not None:
            self.thread.job.cancel()
            self.thread.wait()
        super().closeEvent(event)


def show_dialog():
    global _dialog
    if _dialog is None:
        from pymol.gui import get_qtwindow
        _dialog = AntibodyDialog(get_qtwindow())
    _dialog.show()
    _dialog.raise_()
    _dialog.activateWindow()
    return _dialog
