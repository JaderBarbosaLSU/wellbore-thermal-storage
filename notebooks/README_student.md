# The student copies

These two notebooks are **pinned to the tag `v1.0-2026-10-31`**. Every
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

Change the one line `GITREF = 'v1.0-2026-10-31'` near the top of the
fetch cell to another tag, or to `'notebook-model'` to track the live
line. It repoints every download in the notebook at once.
