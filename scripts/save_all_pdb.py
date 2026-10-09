"""
save_all_pdb.py - Save every loaded PyMOL model as a .pdb file named after its object.

Usage (inside PyMOL):
    run save_all_pdb.py
    save_all_pdb                      # saves into PyMOL's current working directory
    save_all_pdb /path/to/output      # saves into a specific folder (created if missing)
    save_all_pdb ., all_states=0      # save only the current state instead of all states

PyMOL names each object after the file it was loaded from (minus the extension),
so 1abc.cif becomes object "1abc" and is saved as 1abc.pdb.
"""

import os
from pymol import cmd


def save_all_pdb(out_dir=".", all_states=1):
    out_dir = os.path.abspath(os.path.expanduser(str(out_dir)))
    os.makedirs(out_dir, exist_ok=True)
    state = 0 if int(all_states) else -1  # 0 = all states (MODEL/ENDMDL), -1 = current

    objects = cmd.get_object_list("(all)")  # molecular objects only (skips maps, CGOs, etc.)
    if not objects:
        print(" save_all_pdb: no molecular objects loaded.")
        return

    for obj in objects:
        path = os.path.join(out_dir, obj + ".pdb")

        # Warn about things the PDB format can't represent faithfully
        n_atoms = cmd.count_atoms(obj)
        if n_atoms > 99999:
            print(f" Warning: {obj} has {n_atoms} atoms; PDB serial numbers max out at 99999.")
        long_chains = [c for c in cmd.get_chains(obj) if len(c) > 1]
        if long_chains:
            print(f" Warning: {obj} has multi-character chain IDs {long_chains}; "
                  "PDB allows only one character, so they will be truncated.")

        cmd.save(path, obj, state=state, format="pdb")
        print(f" Saved {obj} -> {path}")

    print(f" save_all_pdb: wrote {len(objects)} file(s) to {out_dir}")


cmd.extend("save_all_pdb", save_all_pdb)