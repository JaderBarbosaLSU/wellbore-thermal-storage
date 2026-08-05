# Uploading this to GitHub

I still cannot push from this session — the cloud sandbox's git proxy refuses the
repository, independently of your token. So this has to go up by hand, once.

## Fastest route (one drag, keeps folder structure)

1. Download every file I sent and rebuild this tree on your Desktop:

```
thums-upload/
  pyproject.toml
  thums_core/__init__.py  errors.py  config.py  kpi.py  doe.py  results.py
        stefan.py  well_marched.py  sizing.py  system.py  _legacy_physics.py
        post/__init__.py  effects.py  pareto.py
  tools/extract_legacy.py  freeze_baseline.py  diagnose_closure.py
  notebooks/01_energy_closure.ipynb
  docs/*.md
  verification/reference/thums_v0_baseline.json
```

2. In the repository root, **Add file → Upload files**
3. Drag the *contents* of `thums-upload/` (not the folder itself) into the drop zone.
   GitHub preserves subfolders when you drag folders from your file manager.
4. Commit message: `thums_core v0.2: package skeleton, energy-closure fix, two-constraint sizing`
5. Tag it: https://github.com/JaderBarbosaLSU/wellbore-thermal-storage/releases/new → `v0.2`

Never the web editor — that is what truncated the 2.9 MB notebook to 2 bytes.

## Regenerating `_legacy_physics.py`

It is auto-generated and large. If you would rather not upload it, skip it and run
this once after the rest is up:

```
python tools/extract_legacy.py thums_core/THUMS_Multilayer_BB_param_Jan_12b_converg_improved.ipynb
```

## Then, in Colab

Open the bookmarked link for `notebooks/01_energy_closure.ipynb`:

https://colab.research.google.com/github/JaderBarbosaLSU/wellbore-thermal-storage/blob/main/notebooks/01_energy_closure.ipynb

Runtime → Restart and run all. The first cell installs the package from GitHub, so
after any future change you must **restart**, not just re-run — otherwise Colab keeps
the previously installed version while the notebook text looks current.
