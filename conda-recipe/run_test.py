import os
import sys
import traceback
from subprocess import Popen, check_output

pkg_name = 'chemreac'

# THIS DID NOT WORK (tests run in isolation of source):
# -----------------------------------------------------
# Work around the libm.so bug:
# print(__file__)
# print(check_output(['pwd']))
# print(check_output(['ls', '..']))
#conda_path=check_output(['./scripts/get_conda_path.sh'])
#p = Popen(['bash', './scripts/pruge_conda_shared_obj.sh', './chemreac/_chemreac.so', conda_path, 'libm'])
#p.wait()


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
