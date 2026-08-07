import napari
import numpy as np
import tifffile as tiff

def line_in_time():
    viewer = napari.Viewer()

    path = "/home/clement/Documents/projects/nucleation/3VPCs/251119_#4_30_001_016.vsi - C561.tif"
    img_2dt = tiff.imread(path)

    N = 100
    time_axis = np.arange(50, 50+N, 1)
    p1 = np.arange(10, 10+N, 1)
    p2 = np.arange(100, 100+N, 1)

    p_left = np.stack([time_axis, p1, p1], axis=-1)
    p_right = np.stack([time_axis, p2, p2], axis=-1)
    points = np.stack([p_left, p_right], axis=1)  # (T, 2, 3)

    print(points)

    lines = [
        np.array([
            [163, 10, 10], [163, 100, 100],
        ]),
        np.array([
            [164, 100, 100], [164, 200, 200],
        ]),
    ]

    viewer.add_image(
        img_2dt,
        name="2D+t image"
    )

    viewer.add_shapes(
        points,
        shape_type='line'
    )

    napari.run()

def add_start_coordinate():

    line0 = np.array([
        [163, 10, 10], [163, 100, 100],
    ])
    line1 = np.array([
        [10, 10], [100, 100],
    ])
    start_t = 150
    line = line1
    
    if line.shape[1] == 2:
        line = np.insert(line, 0, start_t, axis=1)
    else:
        line[0, 0] = start_t
        line[1, 0] = start_t

    print(line)


if __name__ == "__main__":
    add_start_coordinate()