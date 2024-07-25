import numpy as np
import matplotlib.pyplot as plt
import imageio


def read_and_plot(filename, imax, jmax):
    """Reads the data from the file and plots it using matplotlib."""
    data = np.fromfile(filename, dtype=np.double).reshape((imax, jmax))
    plt.imshow(data, cmap='hot', interpolation='nearest')
    plt.colorbar()
    plt.title(filename.split('/')[-1])
    return plt.gcf()


def create_gif(imax, jmax, folder_path, output_gif_name, num_files, file_prefix):
    """Creates a gif from the output files."""
    filenames = [f"{folder_path}/{file_prefix}_{i}.dat" for i in range(1, num_files + 1)]
    images = []
    for filename in filenames:
        fig = read_and_plot(filename, imax, jmax)
        fig.savefig('temp.png', bbox_inches='tight')
        images.append(imageio.imread('temp.png'))
        plt.close(fig)

    imageio.mimsave(output_gif_name, images, fps=5)


# Parameters
imax, jmax = 301, 301
num_files = 300
folder_path = 'result_neu'
file_prefix1 = 'output_C1'
file_prefix2 = 'output_C2'
output_name1 = 'animation/neumann_C1.gif'
output_name2 = 'animation/neumann_C2.gif'

# Create the gifs
create_gif(imax, jmax, folder_path, output_name1, num_files, file_prefix1)
create_gif(imax, jmax, folder_path, output_name2, num_files, file_prefix2)
