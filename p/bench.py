#!/usr/bin/env python3
import subprocess
from pathlib import Path
import numpy as np
import sys


s = sys.argv[1]
in_dir = Path("compose/eos/epch")
out_dir = Path(f"output/{s}")

subprocess.run(f"rm -rf {str(out_dir)}", shell=True)
out_dir.mkdir(exist_ok=True, parents=True)

subprocess.run(f"make {s}", shell=True)

N = 21
rtol_min = 1e-8
rtol_max = 1e-3
rtols = np.geomspace(rtol_min, rtol_max, N)

for i, rtol in enumerate(rtols):
    sub_out_dir = out_dir / f"{i}"
    sub_out_dir.mkdir(exist_ok=True, parents=True)
    with open(sub_out_dir / "rtol.txt", "w") as file:
        file.write(str(rtol))
    for in_file in in_dir.glob("*.txt"):
        out_file = sub_out_dir / in_file.name
        subprocess.run([f"./build/{s}", str(in_file), str(out_file), str(rtol)])
        # exit(1)
