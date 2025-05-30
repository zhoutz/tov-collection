import numpy as np


def read_1st_line(filename):
    with open(filename, "r") as file:
        first_line = file.readline().strip().split()
        return np.array([float(x) for x in first_line])

def read_1st_float(filename):
    return read_1st_line(filename)[0]

def cmp_file(gt_fname, cmp_fname):
    ground_truth = np.loadtxt(gt_fname, skiprows=1)
    data = np.loadtxt(cmp_fname, skiprows=1)

    mr_truth = ground_truth
    mr_comp = data

    relative_error = np.abs((mr_truth - mr_comp) / mr_truth)
    mean_relative_error = np.average(relative_error, axis=0)
    max_relative_error = np.max(relative_error, axis=0)
    avgs = read_1st_line(cmp_fname) / data.shape[0]
    avg_time = avgs[0]
    avg_cnt = avgs[1]
    return (mean_relative_error, max_relative_error), avg_time, avg_cnt
