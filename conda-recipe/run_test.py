import os
import sys
import traceback
from subprocess import Popen, check_output

pkg_name = 'chemreac'

print("Importing %s..." % pkg_name)
try:
    import chemreac
except ImportError as e:
    print("Importing %s failed" % pkg_name)
    print(e)
    traceback.print_exc(file=sys.stdout)
    raise e

print("Running tests...")

# Run the full test suite
p = Popen(['py.test', '--pyargs', pkg_name]) # need conftest.py for: '--slow', '--veryslow'
assert p.wait() == os.EX_OK
