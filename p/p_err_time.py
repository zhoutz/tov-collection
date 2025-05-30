from itertools import product
import numpy as np
from ids import mrt_ids
from uts import cmp_file, read_1st_float
import matplotlib.pyplot as plt
from pathlib import Path
import sys


s = sys.argv[1]
src_dir = Path(f"output/{s}")

plt_data_rtol = []
plt_data_cnt = []
plt_data_err = []

for subdir in src_dir.iterdir():
    if not subdir.is_dir() or subdir.name == "0":
        continue

    def f(id):
        gt_fname = src_dir / f"0/{id}.txt"
        cmp_fname = subdir / f"{id}.txt"
        return cmp_file(gt_fname, cmp_fname)

    errs, ts, cnts = zip(*map(f, mrt_ids))
    errs = np.mean(errs, axis=0)
    mean_time = np.mean(ts)
    std_time = np.std(ts)
    mean_cnt = np.mean(cnts)
    std_cnt = np.std(cnts)

    plt_data_rtol.append(read_1st_float(subdir / "rtol.txt"))
    plt_data_err.append(errs)
    plt_data_cnt.append((mean_time, std_time, mean_cnt, std_cnt))

plt_data_err = np.array(plt_data_err)
plt_data_cnt = np.array(plt_data_cnt)

plt.figure(figsize=(6, 7))
plt.subplot(211)
for i, j in product(range(1), range(3)):
    plt.plot(
        plt_data_rtol,
        plt_data_err[:, i, j],
        ["C0", "C1", "C2"][j] + [".", "+"][i],
        label=["", "max"][i] + " " + ["M", "R", "T"][j],
    )
plt.ylabel("mean relative error")
plt.tick_params(
    axis="x",
    labelbottom=False,
)
plt.legend()
plt.xscale("log")
plt.yscale("log")


plt.subplot(212)
plt.plot(plt_data_rtol, plt_data_cnt[:, 0], "g+", label="Mean Time")
plt.ylabel(r"time ($\mu$s)", color="g")
plt.yticks(color="g")
plt.xlabel("rtol")
plt.twinx()
plt.plot(plt_data_rtol, plt_data_cnt[:, 2], "r.", label="Mean Count")
plt.ylabel("# fn eval", color="r")
plt.yticks(color="r")
plt.xscale("log")

plt.suptitle(f"Error and Time Analysis for {s}", fontsize=14)
plt.tight_layout()
plt.subplots_adjust(hspace=0)

if True:
    plt.show()
else:
    saved_fname = f"output/{s}.svg"
    plt.savefig(saved_fname)
    print(f"Saved to ./{saved_fname}")
