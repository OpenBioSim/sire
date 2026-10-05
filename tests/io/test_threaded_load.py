import subprocess
import sys

import sire as sr

_script = """
import sys
import threading

import sire as sr

filename = sys.argv[1]

# Finish importing sire's modules first, since only the loads should overlap.
sr.load(filename, show_warnings=False)

num_threads = 8
barrier = threading.Barrier(num_threads)
errors = []


def load():
    barrier.wait()
    try:
        for _ in range(5):
            mols = sr.load(filename, show_warnings=False)
            mols.molecules("not water")
    except Exception as e:
        errors.append(e)


threads = [threading.Thread(target=load) for _ in range(num_threads)]
for thread in threads:
    thread.start()
for thread in threads:
    thread.join()

if errors:
    raise errors[0]
"""


def test_threaded_load(tmpdir, kigaki_mols):
    # Run in a separate process, since the failure is a crash.
    d = tmpdir.mkdir("test_threaded_load")
    f = sr.save(kigaki_mols, d.join("system"), format=["PRM7"])

    proc = subprocess.run(
        [sys.executable, "-c", _script, f[0]], capture_output=True, text=True
    )

    assert proc.returncode == 0, proc.stderr
