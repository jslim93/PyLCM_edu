import pytest
from PyLCM.mixing import LEMMixing


def test_lem_backend_is_a_clear_stub():
    lem = LEMMixing()
    with pytest.raises(NotImplementedError, match="Phase 3b"):
        lem.apply([], 288.0, 0.009, 90000.0, 500.0, 1.0, 1.0, 100.0)
