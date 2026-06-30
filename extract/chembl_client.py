"""ChEMBL REST API client: session management, retry/backoff, pagination.

Implemented in milestone M1 (feat/extract-raw). Kept import-light so the
package imports cleanly without the ``extract`` extra installed (e.g. in M0 CI).
"""


def get_client():
    """Return a configured ChEMBL web resource client."""
    raise NotImplementedError("Implemented in M1 (feat/extract-raw).")
