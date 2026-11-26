'''
The menu of the toolbar's Plugins button: every entry of the Plugin
menu in one list, under a heading for what the tools are for.

The entries are taken from the Plugin menu each time the list is shown,
so a plugin installed later (Plugin Manager) or a legacy one is listed
too. A plugin no group names is listed under "Miscellaneous": name a
new plugin's menu label in GROUPS to give it a place.
'''

from pymol.Qt import QtWidgets

# (heading, the plugins' labels in the Plugin menu, in the order shown)
GROUPS = [
    ('Sequence Tools', [
        'Sequence Comparison',
        'Modify Sequences',
    ]),
    ('Structure Tools', [
        'Structure Comparison',
        'Structure Prediction',
        'Structure Preparation',
        'Interface Analysis (PISA-like)',
        'Alignment/Superposition',
    ]),
    ('Clustering', [
        'RMSD Clustering',
        'Interface RMSD Clustering',
    ]),
    ('Antibody Tools', [
        'Antibody Builder',
        'CDR Loops',
        'Epitope Prediction',
        'Paratope Prediction',
    ]),
    ('Ligand Tools', [
        'Binding Pocket',
    ]),
    ('Physicochemical Tools', [
        'Aggregation Prediction',
        'APBS Electrostatics',
    ]),
]
MISCELLANEOUS = 'Miscellaneous'

# not tools: listed last, under no heading
LAST = ['Plugin Manager', 'Initialize Plugin System']


def _key(label):
    return label.replace('&', '').strip().lower()


def grouped(labels):
    '''
    Sort menu labels into the groups.

    :param labels: labels of the Plugin menu, in its order
    :return: ([(heading, [label])] without empty groups, [label] of LAST)
    '''
    labels = list(dict.fromkeys(labels))
    known = {_key(label): label for label in labels}
    last = [known.pop(_key(name)) for name in LAST if _key(name) in known]
    groups = []
    for heading, names in GROUPS:
        members = [known.pop(_key(name)) for name in names
                   if _key(name) in known]
        if members:
            groups.append((heading, members))
    # what is left, as the Plugin menu orders it
    others = [label for label in labels if _key(label) in known]
    if others:
        groups.append((MISCELLANEOUS, others))
    return groups, last


def source_menu(window):
    '''
    The main window's Plugin menu (None without one)
    '''
    menus = getattr(window, 'menudict', None) or {}
    # "PluginQt" once the plugins are loaded: "Plugin" is its Legacy
    # Plugins submenu from then on
    return menus.get('PluginQt') or menus.get('Plugin')


def plugin_actions(source):
    '''
    The entries of a menu and of its submenus (Legacy Plugins), without
    the separators.

    :return: [QAction]
    '''
    found = []
    for action in (source.actions() if source is not None else ()):
        if action.isSeparator():
            continue
        submenu = action.menu()
        if submenu is not None:
            found.extend(plugin_actions(submenu))
        elif action.text():
            found.append(action)
    return found


def fill(menu, source):
    '''
    List the entries of the Plugin menu "source" in "menu", grouped.
    The actions are the Plugin menu's own, so an entry does the same
    in either place.

    :return: [(heading, [label])] as listed
    '''
    for action in menu.actions():
        menu.removeAction(action)
        if action.property('plugin_heading') or action.isSeparator():
            action.deleteLater()
    actions = {}
    for action in plugin_actions(source):
        actions.setdefault(action.text(), action)
    groups, last = grouped(actions)
    for heading, members in groups:
        if menu.actions():
            menu.addSeparator()
        title = menu.addAction(heading)
        title.setEnabled(False)
        title.setProperty('plugin_heading', True)
        for label in members:
            menu.addAction(actions[label])
    if last and menu.actions():
        menu.addSeparator()
    for label in last:
        menu.addAction(actions[label])
    if not menu.actions():
        menu.addAction('No plugins').setEnabled(False)
    return groups


def create(window, parent=None):
    '''
    The menu, filled anew from the window's Plugin menu when shown.
    '''
    menu = QtWidgets.QMenu(parent)
    menu.setObjectName('plugins_menu')

    def update():
        fill(menu, source_menu(window))

    menu.aboutToShow.connect(update)
    return menu
