# Based upon:
# http://stackoverflow.com/a/18143318/790973

"""
Usage:

.. exec::
   echo "Shell code!"

"""

import codecs
from os.path import basename
import subprocess
import sys
import traceback

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

from sphinx.util.compat import Directive
from docutils import nodes, statemachine


class ExecDirective(Directive):
    # See: http://docutils.sourceforge.net/docs/howto/rst-directives.html
    """Execute the specified command and insert the output into the document"""
    has_content = True
    required_arguments = 0
    optional_arguments = 0
    final_argument_whitespace = False
    option_spec = {}

    def run(self):
        tab_width = self.options.get(
            'tab-width', self.state.document.settings.tab_width)
        source = self.state_machine.input_lines.source(
            self.lineno - self.state_machine.input_offset - 1)
        try:
            output = []
            for cmd in self.content:
                output.append(subprocess.check_output(cmd, shell=True))
            lines = statemachine.string2lines(
                '\n'.join([codecs.decode(elem, 'utf-8') for elem in output]),
                tab_width, convert_whitespace=False)
            self.state_machine.insert_input(lines, source)
            return []
        except Exception:
            exc = sys.exc_info()
            return [nodes.error(None, nodes.paragraph(
                text="Unable to run command at %s:%d:" % (
                    basename(source), self.lineno)
            ), nodes.paragraph(text=str(exc[1]) + '\n'.join(
                [] + traceback.format_tb(exc[2]))))]


def setup(app):
    app.add_directive('exec', ExecDirective)
