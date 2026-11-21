"""The package as a user starts it: what `import muoscene` does to the process."""

import os
import subprocess
import sys


def test_mulder_finds_its_cache_without_home():
    """Windows sets no `HOME`, and Mulder 0.3.7 resolves its cache through it alone.

    The first run of the tests on Windows, 21 November 2026: 144 tests failed with "bad cache
    (could not resolve $HOME)". A fresh interpreter without `HOME`: after `import muoscene` the
    variable is the home directory of the user, Mulder has a cache directory under it and prints
    no panic. Nothing is compiled here, so nothing is written there.
    """
    code = (
        "import os, pathlib, muoscene, mulder\n"
        "home = pathlib.Path.home()\n"
        "assert os.environ['HOME'] == str(home), os.environ.get('HOME')\n"
        "assert mulder.config.CACHE == home / '.cache' / 'mulder', mulder.config.CACHE\n"
    )
    environment = {name: value for name, value in os.environ.items() if name != "HOME"}
    run = subprocess.run(
        [sys.executable, "-c", code], env=environment, capture_output=True, text=True, timeout=120
    )
    assert run.returncode == 0, run.stderr
    assert "panicked" not in run.stderr


def test_a_home_that_is_set_is_left_alone():
    code = "import os, muoscene; print(os.environ['HOME'])"
    environment = dict(os.environ, HOME=os.sep.join(["", "somewhere", "else"]))
    run = subprocess.run(
        [sys.executable, "-c", code], env=environment, capture_output=True, text=True, timeout=120
    )
    assert run.returncode == 0, run.stderr
    assert run.stdout.strip() == environment["HOME"]
