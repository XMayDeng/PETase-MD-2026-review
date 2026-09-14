"""Run inside ChimeraX to project exact residue anchors before a PNG save.

No coordinates or cameras are modified. The final compositor draws crisp
screen-plane labels; the molecular render itself contains no label textures.
"""

import json
import os
from pathlib import Path

import numpy as np
from chimerax.core.commands import run
from chimerax.label.label3d import ObjectLabels


WORK = Path(os.environ["PETASE_RENDER_WORK"])
NAMES = {1: "IsPETase", 2: "TfCut1", 3: "LCC", 4: "FoCut5a", 5: "HiC", 6: "PHL7"}
ANCHOR_ATOMS = {
    "SER": {"OG"}, "ASP": {"CG", "OD1", "OD2"},
    "HIS": {"CG", "ND1", "CD2", "CE1", "NE2"},
    "TRP": {"CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"},
}
view = session.main_view
camera = view.camera
image_size = (2100, 1800)
near_far = view.near_far_distances(camera, None)
projection = camera.projection_matrix(near_far, None, image_size)
inverse_camera = camera.position.inverse().opengl_matrix()


def screen_fraction(scene_coordinates):
    point = np.append(scene_coordinates, 1.0) @ inverse_camera @ projection
    point = point[:3] / point[3]
    return [float((point[0] + 1) / 2), float((point[1] + 1) / 2)]


records = []
for model in session.models.list(type=ObjectLabels):
    for label in model.labels():
        residue = label.object
        atoms = [atom for atom in residue.atoms if atom.name in ANCHOR_ATOMS[residue.name]]
        if not atoms:
            raise ValueError(f"No anchor atoms for {residue}")
        anchor = np.mean([atom.scene_coord for atom in atoms], axis=0)
        selected_atoms = [atom for atom in residue.atoms if atom.element.number > 1]
        records.append({
            "label": label.text, "model_id": residue.structure.id[0],
            "chain_id": residue.chain_id, "pdb_residue_number": int(residue.number),
            "residue_name": residue.name,
            "anchor_atom_names": [atom.name for atom in atoms],
            "anchor_scene_angstrom": anchor.tolist(),
            "anchor_screen_fraction": screen_fraction(anchor),
            "heavy_atom_screen_fractions": [screen_fraction(atom.scene_coord) for atom in selected_atoms],
        })

model_ids = {record["model_id"] for record in records}
if len(model_ids) != 1:
    raise ValueError(f"Expected one currently labelled structure, found {model_ids}")
name = NAMES[model_ids.pop()]
report = {
    "protein": name, "image_size_px": list(image_size),
    "camera_name": camera.name, "camera_position": camera.position.matrix.tolist(),
    "camera_projection_matrix": np.asarray(projection).tolist(),
    "model_scene_position": residue.structure.scene_position.matrix.tolist(),
    "records": records,
}
(WORK / f"Figure_01_{name}_label_projection.json").write_text(json.dumps(report, indent=2) + "\n")
run(session, "label delete")
assert np.array_equal(camera.position.matrix, np.asarray(report["camera_position"]))
assert np.array_equal(residue.structure.scene_position.matrix, np.asarray(report["model_scene_position"]))
