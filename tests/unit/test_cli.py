import pytest

from muoscene import __version__
from muoscene.cli import COMMANDS, main


def test_no_command_prints_help(capsys):
    assert main([]) == 0
    out = capsys.readouterr().out
    for command in COMMANDS:
        assert command in out


def test_version(capsys):
    with pytest.raises(SystemExit) as stop:
        main(["--version"])
    assert stop.value.code == 0
    assert __version__ in capsys.readouterr().out


@pytest.mark.parametrize("command", sorted(COMMANDS))
def test_commands_parse(command, capsys):
    assert main([command, "run.yaml"]) == 2
    assert command in capsys.readouterr().err


def test_unknown_command():
    with pytest.raises(SystemExit) as stop:
        main(["no-such-command"])
    assert stop.value.code == 2
