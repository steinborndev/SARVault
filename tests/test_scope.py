"""Tests for the global scope cross-filter (dashboard/scope.py).

Covers the pure target-reconcile helper and, via Streamlit's AppTest harness, the
stickiness of the scope controls: adding or removing a target must survive the rerun
the edit triggers (the bug these tests guard against was edits silently snapping back).
"""

import pytest

from dashboard import scope


def test_reconcile_targets_drops_stale_and_keeps_order():
    options = ["Tubulin", "DNA topoisomerase 1", "HDAC1"]
    # a stale entry is dropped, valid ones kept in their given order
    assert scope.reconcile_targets(["HDAC1", "GHOST", "Tubulin"], options) == ["HDAC1", "Tubulin"]
    # empty / fully-stale selections collapse to every option (never filter to nothing)
    assert scope.reconcile_targets([], options) == options
    assert scope.reconcile_targets(["GHOST"], options) == options
    assert scope.reconcile_targets(None, options) == options


def _scope_harness():
    from dashboard import scope

    targets = ["Tubulin", "DNA topoisomerase 1", "DNA topoisomerase 2-alpha", "HDAC1"]
    scope.render(targets)


def _targets(at):
    return at.multiselect[0].value


def test_scope_controls_persist_edits():
    AppTest = pytest.importorskip("streamlit.testing.v1").AppTest

    at = AppTest.from_function(_scope_harness).run()
    assert not at.exception
    assert "DNA topoisomerase 1" in _targets(at)

    # remove a target -> the edit must stick across the rerun it triggers
    at.multiselect[0].unselect("DNA topoisomerase 1").run()
    assert "DNA topoisomerase 1" not in _targets(at)

    # an unrelated rerun must NOT re-seed the widget back to all targets
    at.run()
    assert "DNA topoisomerase 1" not in _targets(at)

    # add it back -> also sticks
    at.multiselect[0].select("DNA topoisomerase 1").run()
    assert "DNA topoisomerase 1" in _targets(at)

    # a target dropped from the warehouse is reconciled out without crashing
    at.session_state["scope_targets"] = ["Tubulin", "GHOST_TARGET"]
    at.run()
    assert not at.exception
    assert _targets(at) == ["Tubulin"]
