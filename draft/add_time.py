import numpy as np

n_times = 447
n_points = 29
time_start = 17

base_table = np.zeros((n_times, n_points, 2), dtype=np.float32)

time_range = np.arange(time_start, time_start + n_times)
time_stack = np.repeat(time_range[:, np.newaxis], n_points, axis=1)
time_stack = time_stack[..., np.newaxis]
print(time_stack.shape)