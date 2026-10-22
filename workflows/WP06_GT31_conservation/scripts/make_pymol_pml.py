"""Create a self-contained PyMOL command file with embedded class PDBs.

The generated PML does not depend on PyMOL's working directory. It restores
core coloring, not the historical superposition or camera orientation.
"""
from pathlib import Path


def build(reference_pdb: Path, query_pdb: Path, output: Path) -> None:
    lines = ['# Generated core-coloring view; PDB data are embedded.',
             '# B factors encode classes, not pLDDT or experimental displacement.',
             'python', 'from pymol import cmd']
    for name, path in [('Q18515_core', reference_pdb), ('F2Y1U4_core', query_pdb)]:
        lines.extend([f'cmd.delete({name!r})', f'cmd.read_pdbstr({path.read_text()!r}, {name!r})'])
    lines.extend([
        'python end',
        'hide everything, Q18515_core or F2Y1U4_core',
        'show cartoon, Q18515_core or F2Y1U4_core',
        'set cartoon_color, -1, Q18515_core or F2Y1U4_core',
        'color grey80, Q18515_core or F2Y1U4_core',
        'set_color gt31_strict, [0.000, 0.486, 0.569]',
        'set_color gt31_supported, [0.651, 0.808, 0.890]',
        'select gt31_strict_positions, (Q18515_core or F2Y1U4_core) and b > 99',
        'select gt31_supported_positions, (Q18515_core or F2Y1U4_core) and b > 59 and b < 61',
        'color gt31_strict, gt31_strict_positions',
        'color gt31_supported, gt31_supported_positions',
        'bg_color white', 'set ray_opaque_background, off', 'set antialias, 2',
        'deselect', 'zoom Q18515_core or F2Y1U4_core',
        '# No structural superposition or manuscript camera view is imposed.',
    ])
    output.write_text('\n'.join(lines) + '\n')
