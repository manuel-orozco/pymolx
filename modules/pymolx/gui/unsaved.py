'''
Before quitting: offer to save the session if it has unsaved changes
(pymolx.unsaved), with Save, Don't Save and Cancel
'''
import os

from pymol.Qt import QtWidgets

from .. import branding, unsaved

QMessageBox = QtWidgets.QMessageBox


def question(window, _self):
    '''
    The question, not shown yet (tests answer it themselves)
    '''
    box = QMessageBox(window)
    box.setObjectName('unsaved_question')
    box.setIcon(QMessageBox.Icon.Question)
    box.setWindowTitle('Quit %s' % branding.NAME)
    box.setText('Save the current session before quitting?')
    session = _self.get('session_file')
    box.setInformativeText(
        'Changes made since %s was saved will be lost if you don\'t save '
        'them.' % os.path.basename(session) if session else
        'The session has not been saved; what it holds will be lost if you '
        'don\'t save it.')
    box.setStandardButtons(QMessageBox.StandardButton.Save |
                           QMessageBox.StandardButton.Discard |
                           QMessageBox.StandardButton.Cancel)
    box.button(QMessageBox.StandardButton.Discard).setText('Don\'t Save')
    box.setDefaultButton(QMessageBox.StandardButton.Save)
    box.setEscapeButton(QMessageBox.StandardButton.Cancel)
    return box


def answer(window, _self, button):
    '''
    Act on the answer: True to quit (saved, or not wanted), False to stay
    (cancelled, or the session wasn't saved after all)
    '''
    if button == QMessageBox.StandardButton.Discard:
        return True
    if button == QMessageBox.StandardButton.Save:
        # to the session file it came from or was saved as, else ask where
        window.session_save()
        return not unsaved.has_unsaved(_self)
    return False


def confirm(window, _self):
    '''
    True if quitting may go on: nothing unsaved, or the user saved the
    session or chose not to
    '''
    if not unsaved.has_unsaved(_self):
        return True
    box = question(window, _self)
    box.exec()
    return answer(window, _self, box.standardButton(box.clickedButton()))
