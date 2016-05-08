#  import matplotlib.pyplot as plt
#  import seaborn as sns

import matplotlib.pyplot as plt
from os.path import join
# from helpers.config import Config


class BasePlotter:
    """
    Abstract class. Use PCAPlotter and the like, not this one.
    """
    def __init__(self):
        # FIXME make this work
        # self.colors = Config('colors')
        pass

    def savefig(self, filename):
        if not filename.endswith('.png'):
            filename += '.png'
        filepath = join(self.base_dir, filename)
        plt.savefig(filepath, dpi=300, bbox_inches='tight')
        print('Saved at -> ' + filepath)
        return filepath
