import pytest
from PyLCM.entrainment import basic_entrainment


def test_entrainment_blocked_by_default():
    # The experimental gate must raise before any physics runs.
    with pytest.raises(RuntimeError, match="EXPERIMENTAL"):
        basic_entrainment(1.0, 100.0, 285.0, 0.008, 95000.0, 0.01, None, None)


def test_entrainment_opt_in_passes_the_gate():
    # With experimental=True the gate is passed; it then fails LATER on the
    # undefined env profile (a separate known limitation), not on the gate.
    with pytest.raises(Exception) as exc:
        basic_entrainment(1.0, 100.0, 285.0, 0.008, 95000.0, 0.01, None, None,
                          experimental=True)
    assert "EXPERIMENTAL" not in str(exc.value)
