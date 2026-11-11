# The student copies

These two notebooks are **pinned to the tag `v1.1-2026-11-11`**. Every
file they download — the slot coordinates, the run cache, the cached
field solve — comes from that tag, which never moves. Nothing pushed to
`notebook-model` afterwards can change what they do.

That is the point: you can experiment here without any risk of being
moved under, and without any risk of disturbing the working line.

## What you cannot break

- Opening a notebook from GitHub in Colab gives you a **read-only view**.
  You cannot save back to this repository from it.
- `File → Save a copy in Drive` makes **your own** copy in your Google
  Drive. Edit that one freely; it is yours and it touches nothing here.
- The `notebook-model` branch is the live line and is not affected by
  anything you do on `student`.

## If you want your work kept

Commit it to this `student` branch, or open a pull request against it.
Do not push to `notebook-model`.

## To move to a newer snapshot

Change the one line `GITREF = 'v1.1-2026-11-11'` near the top of the
fetch cell to another tag, or to `'notebook-model'` to track the live
line. It repoints every download in the notebook at once.

**`BOREAS_model.ipynb` sets `GITREF` in two cells (4 and 6) — change
both, or `tools/student_check.py` will tell you they disagree.**

## What changed from v1.1-2026-11-11

v1.0 had two defects, and both are fixed here:

- `BOREAS_model.ipynb` carried a stray line in the pack-fetch cell and
  **raised `SyntaxError: invalid syntax`**;
- only one of its two `GITREF` assignments was pinned to the tag, so the
  snapshot was advertised as frozen while fetching part of its inputs
  from the live branch.

`tools/student_check.py` now checks both before any snapshot is cut:

    python3 tools/student_check.py notebooks/BOREAS_model.ipynb \
        notebooks/BOREAS_parametric.ipynb --expect v1.1-2026-11-11

The model itself also advances to pack v2.8. No published number moves,
but the water-loop state points are **renumbered**: 1 is now the field
outlet and 4 the field inlet, on both half-cycles.
