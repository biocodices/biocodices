import pandas as pd
import matplotlib.pyplot as plt

from biocodices.plotters.association_test_plotter import AssociationTestPlotter


class AssociationTest:
    def __init__(self, dataset):
        self.dataset = dataset
        self.plotter = AssociationTestPlotter(self, self.dataset.dir)

    def __repr__(self):
        return '<{} for {}>'.format(self.__class__.__name__,
                                    self.dataset.label)

    def run(self):
        """
        Run a battery of association tests per SNP and store the output in
        .result. FIXME: should be improved.
        """
        fn1 = self.dataset.plink.assoc(adjust=False)
        df1 = pd.read_table(fn1, sep='\s+').set_index('SNP')
        fn2 = self.dataset.plink.assoc(adjust=True)
        df2 = pd.read_table(fn2, sep='\s+').set_index('SNP')

        df = df1.join(df2, rsuffix='_dupe')
        fn1 = self.dataset.plink.assoc(adjust=False)
        df1 = pd.read_table(fn1, sep='\s+').set_index('SNP')
        fn2 = self.dataset.plink.assoc(adjust=True)
        df2 = pd.read_table(fn2, sep='\s+').set_index('SNP')

        df = df1.join(df2, rsuffix='_dupe')
        df = df.drop([col for col in df.columns if '_dupe' in col], axis=1)
        self.result = df

        fn3 = self.dataset.plink.model()
        df3 = pd.read_table(fn3, sep='\s+').set_index('SNP')
        self.available_models = []
        for model, df in df3.groupby('TEST'):
            model_name = model.lower()
            self.available_models.append(model_name)
            setattr(self, model_name, df)
            self.result[model] = df['P']

        self.available_tests = self.result.ix[:, 'UNADJ':].columns

    def plot(self, ax=None, tests_to_plot=None):

        if ax is None:
            _, ax = plt.subplots(figsize=(15, 5))

        if tests_to_plot is None:
            tests_to_plot = self.available_tests

        self.plotter.draw_ax(ax, tests_to_plot)
        return ax

    def plot_models(self, ax=None):
        if ax is None:
            _, ax = plt.subplots(figsize=(15, 5))

        for model in self.available_models:
            df = getattr(self, model)
            df.rename(columns={'P': model})
            self.plotter.draw_ax(ax, tests_to_plot=[model])

        return ax
