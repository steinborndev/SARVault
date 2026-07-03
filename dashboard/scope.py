"""Global scope cross-filter: stateful Streamlit controls + a pure reconcile helper.

The scope popover in the header drives every page. Its controls must be *sticky*:
adding or removing a target has to survive the rerun the edit itself triggers. The
earlier implementation derived each widget's ``default`` from the previously stored
selection, which churned a keyless widget's identity on that rerun and let Streamlit
re-seed it from the stale default — so an edit didn't always take. Here each control
owns a stable session-state key and takes no ``default``, so the widget owns its value
across reruns and every edit persists.
"""

import streamlit as st

_APPROVAL = ["all", "approved", "research"]


def reconcile_targets(stored, options):
    """Keep only stored targets the warehouse still offers; fall back to all options.

    Pure (no Streamlit) so it is unit-testable. An empty or fully-stale selection
    collapses to every option, so the scope never silently filters everything out and
    the multiselect — whose value is set from this — never receives an option it does
    not know (which would raise).
    """
    options = list(options)
    valid = [t for t in (stored or []) if t in options]
    return valid or options


def render(target_names):
    """Render the scope popover and return (and store) the resolved scope dict."""
    prev = st.session_state.get("scope", {})
    st.session_state.setdefault("scope_targets", prev.get("targets", target_names))
    st.session_state.setdefault("scope_approval", prev.get("approval", "all"))
    st.session_state.setdefault("scope_min_pchembl", float(prev.get("min_pchembl", 0.0)))
    st.session_state.setdefault("scope_structure_only", bool(prev.get("structure_only", False)))
    st.session_state["scope_targets"] = reconcile_targets(
        st.session_state["scope_targets"], target_names
    )

    _, scope_col = st.columns([11, 1])
    with scope_col.popover("Scope", width="stretch"):
        st.multiselect("Targets", target_names, key="scope_targets")
        st.radio("Approval", _APPROVAL, horizontal=True, key="scope_approval")
        st.slider("Min best pChEMBL", 0.0, 12.0, step=0.5, key="scope_min_pchembl")
        st.toggle("Only with 3D crystal structure", key="scope_structure_only")

    scope = {
        "targets": st.session_state["scope_targets"] or list(target_names),
        "approval": st.session_state["scope_approval"],
        "min_pchembl": st.session_state["scope_min_pchembl"],
        "structure_only": st.session_state["scope_structure_only"],
    }
    st.session_state["scope"] = scope
    return scope
