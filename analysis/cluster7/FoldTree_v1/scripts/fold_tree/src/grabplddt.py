import json
import AFDB_tools
import glob
import numpy as np

class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)

path = ''.join( [ p + '/' for p in snakemake.input[0].split('/')[:-1]]) +'structs/*.pdb'
all_structs ={}
for struct in glob.glob(path):
    all_structs[struct] = AFDB_tools.descr(struct)
with open(snakemake.output[0], 'w' ) as rfout:
    rfout.write(json.dumps(all_structs, cls=NpEncoder))
