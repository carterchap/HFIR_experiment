# gdml/ — detailed chamber geometry

Auto-converted (mostly) CAD geometry for the JHU DAMIC-M box chamber, loaded by
[`DetectorConstruction.cc`](../src/DetectorConstruction.cc). Each `.gdml` file
here has its own header comment (right after the `<gdml ...>` opening tag)
documenting its source, conversion mode, and named volumes — open a file and
scroll to the top before digging into the vertex data below it.

## Layout

```
gdml/
  chamber/     the chamber shell/door
  cryostat/    cryocooler + cold tip/strap/copper block
  vacuum_spacers/
  storage_boxes/
  ccd/         hand-written CCD sensor model
```

| File | Physical part | Volumes | Source | Mode |
|---|---|---|---|---|
| `chamber/box_chamber_jhu_detailed.gdml` | Chamber shell, 5 port flanges, mounting bosses, misc door hardware | 19 | `Box_Chamber_Config_16937.stp.step` | hierarchy |
| `chamber/thedoor.gdml` | Access door: 4 hinges, 12-bolt cover plate, latch/handle | 1 (merged mesh) | Shapr3D STL (`THE_DOOR.zip`) | STL, single mesh |
| `chamber/vacuum_space.gdml` | Interior vacuum cavity (incl. 6 port necks) — gets the Vacuum material, acts as mother volume for parts inside it | 1 (merged mesh) | `Vacuum_Space.step` (loose file, not in the zip) | flat |
| `cryostat/cyro_cooler.gdml` | Cryocooler: pressure vessel, air fins, cold-tip connector | 6 | `Cyro Cooler.step` | hierarchy |
| `cryostat/cold_tip.gdml` | Cold finger | 1 (merged mesh) | `Cold Tip.step` | flat |
| `cryostat/copper_block_detailed.gdml` | Thermal mass block | 1 (merged mesh) | `Copper Block.step` | flat |
| `cryostat/jhu1_cold_strap.gdml` | Cold strap (cold tip → copper block) | 1 (merged mesh) | `JHU1_cold_strap.step` | flat |
| `vacuum_spacers/vacuum_spacer_0.gdml` | Vacuum spacer + DSub feedthrough | 3 | `vacuum_spacer_Rev-.stp.step` | hierarchy |
| `vacuum_spacers/vacuum_spacer_1.gdml` | Vacuum spacer | 2 | `vacuum_spacer_Rev-.stp (1).step` | hierarchy |
| `vacuum_spacers/vacuum_spacer_2.gdml` | Vacuum spacer | 2 | `vacuum_spacer_Rev-.stp (2).step` | hierarchy |
| `storage_boxes/storage_box_1.gdml` | CCD storage box 1/3 | 7 | `6K x 1k Storage Box Assembly...(1).step` | hierarchy |
| `storage_boxes/storage_box_2.gdml` | CCD storage box 2/3 | 6 | `...(2).step` | hierarchy |
| `storage_boxes/storage_box_3.gdml` | CCD storage box 3/3 | 6 | `...(3).step` | hierarchy |
| `ccd/ccd_module.gdml` | 6k×1k CCD sensor (11 layers) | 11 | hand-written | n/a |

"Volumes" counts each file's real named parts (not counting the throwaway
container box every conversion wraps them in — see `UnwrapAndPlaceGDMLContainer`
in `DetectorConstruction.cc`, which discards that box on load).

## Placement model

Every part's absolute chamber-frame position is baked directly into its own
tessellated vertex coordinates (not into a placement transform — every
`<physvol>` in these files uses an identity rotation/translation).
`DetectorConstruction.cc` places each file's top volume at the world origin
with no rotation, which is correct *because* the coordinates are already
absolute. This was verified for every file by comparing bounding boxes against
the chamber shell (vacuum spacers land exactly on port flanges, storage boxes
land inside the shell, the door sits flush on the shell face, etc.) — if you
add a new part and it doesn't line up, suspect the STEP/STL export's
coordinate frame before suspecting the placement code.

## Vacuum nesting

The World volume is `G4_AIR`. Parts confirmed by bounding-box containment
to sit entirely inside `chamber/vacuum_space.gdml`'s cavity are placed as
*real Geant4 daughters* of that volume (given the Vacuum material by
`LoadGDMLSolidWithMaterial` in `DetectorConstruction.cc`), not just
co-located with it — so the empty space around them is correctly Vacuum,
not Air:

- `cryostat/copper_block_detailed.gdml`
- `cryostat/jhu1_cold_strap.gdml`
- `storage_boxes/storage_box_{1,2,3}.gdml` (and the CCDs placed inside them)

Everything else that touches the chamber interior (`cryostat/cyro_cooler.gdml`,
`cryostat/cold_tip.gdml`, `chamber/thedoor.gdml`, the 3 `vacuum_spacers/`) sits
at least partially *outside* the modeled cavity — `cold_tip` in particular
pokes ~32mm beyond it in Z — so they stay as direct World daughters. Geant4
fatally panics ("daughter entirely outside mother") if you nest something
that isn't fully contained, so if you add a new part here, check its bbox
against `vacuum_space.gdml`'s before nesting it, don't assume.

## Editing a part

- **Hierarchy-mode files** (most of them): each real component is its own
  named `<volume>`/`<physvol>` pair — find it by name (see the file's header)
  and edit its `<solid>`, swap its `materialref`, etc.
- **Flat-mode / STL files** (`cold_tip`, `copper_block_detailed`,
  `jhu1_cold_strap`, `thedoor`): the whole part is one merged
  `<tessellated>` solid — there's no sub-part boundary to edit, only the
  whole mesh.
- Every file's top-level `world_lv`/`world_lv_flat`-style container was
  renamed to `<Label>_world_lv[_flat]` on purpose (see the header comments) —
  don't rename it back to plain `world_lv`, or it'll collide with another
  file's top volume in Geant4's global `G4LogicalVolumeStore` and the wrong
  geometry will silently get placed (this happened once already; see git
  history / `DetectorConstruction.cc`'s comments on
  `LoadAndPlaceGDMLAssembly`).

## Regenerating from CAD source

Don't hand-edit the auto-generated files if you're just going to reconvert
from CAD anyway — use `../scripts/regenerate_gdml.py`:

```bash
python3 scripts/regenerate_gdml.py <STEP_ZIP> [DOOR_STL_ZIP] [VACUUM_SPACE_STEP]
```

or via CMake: `cmake --build build --target regenerate_gdml -- STEP_ZIP=... DOOR_ZIP=... VACUUM_ZIP=...`

See that script's module docstring for the full per-file conversion notes
(why each file uses hierarchy vs. flat mode, and why the door specifically
needs a separate STL source rather than the STEP zip).

`ccd/ccd_module.gdml` is hand-written and is **not** touched by
`regenerate_gdml.py` — edit it directly.
