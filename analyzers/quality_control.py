import pandas as pd
import matplotlib.pyplot as plt

from plotters.quality_control_plotter import QualityControlPlotter


class QualityControl:
    def __init__(self, dataset):
        self.dataset = dataset
        self.plotter = QualityControlPlotter(self, self.dataset.dir)

    def __repr__(self):
        return '<{} for {}>'.format(self.__class__.__name__,
                                    self.dataset.label)

    def run(self):
        self.result = {}

        tests = {}
        tests['missing'] = self.dataset.plink.test_missing()

        for test_name, test_filename in tests.items():
            df = pd.read_table(test_filename, sep='\s+').set_index('SNP')
            self.result[test_name] = df

    def plot(self):
        pass
