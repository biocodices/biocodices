import pandas as pd
import subprocess

from os.path import expanduser, join, isfile, getsize
from shutil import copyfile

from .base_pca import BasePCA
from biocodices.helpers.language import percentage_fmt
from biocodices.helpers.config import Config


class SmartPCA(BasePCA):
    _EXECUTABLE = Config('executables')['smartpca']

    def __init__(self, dataset):
        self.dataset = dataset
        self._evecfile = self._output_filepath('pca.evec')
        self._evalfile = self._output_filepath('pca.eval')
        self._logfile = self._output_filepath('pca.log')

    def run(self, overwrite=False, args={}):
        analysis_files_dont_exist = (not isfile(self._evecfile)
                                     or not isfile(self._evalfile))

        if analysis_files_dont_exist or overwrite:
            self._call_smartpca(args)
        else:
            if getsize(self._evecfile) == 0 or getsize(self._evalfile) == 0:
                self._call_smartpca(args)

        self._read_the_results_files()

    def _call_smartpca(self, args={}):
        self.args = {**self.arguments(), **args}
        parfile_path = self._create_parameters_file()
        command = '{} -p {}'.format(self._EXECUTABLE, parfile_path)
        with open(self._output_filepath('pca.log'), 'w+') as logfile:
            subprocess.call(command.split(' '), stdout=logfile)

    def _read_the_results_files(self):
        result = pd.read_table(self._evecfile, sep="\s+", header=None,
                               skiprows=1)
        self.result = self._parse_evec_file(result)
        self.explained_variance = self._read_eval_file(self._evalfile)
        self.write_result_csvs()  # Useful for external use, e.g. d3

        # TODO: read the interesting info in the log file
        # self.extra_info = self._read_log()

    def arguments(self):
        # See ./POPGEN/README in the eigensoft package for a description
        # about each of these parameters and some extra ones.
        args = {
            'genotypename': self.dataset.ped,
            'snpname': self._create_pedsnp(),
            'indivname': self._create_pedind(),
            'numoutevec': 15,  # PCs to take
            'evecoutname': self._evecfile,
            'evaloutname': self._evalfile,
            'altnormstyle': 'NO',
            'numoutlieriter': 5,  # max outlier removal iterations
            'numoutlierevec': 10,  # PCs along which to remove outliers
            'outliersigmathresh': 6,  # min standard deviations of outliers
            'missingmode': 'NO',  # set to YES for 'informative missingness'
            'fastmode': 'NO',
            'outliermode': 1,
        }
        return args

    def _create_parameters_file(self):
        parfile_path = join(self.dataset.label_path + '.pca.par')
        with open(parfile_path, 'w+') as parfile:
            for argname, argvalue in self.args.items():
                parfile.write('{}: {}\n'.format(argname, argvalue))
        return parfile_path

    def _create_pedsnp(self):
        # .pedsnp format is exactly the same as .bim, but smartpca needs
        # the file to have that extension.
        pedsnp_filepath = self.dataset.label_path + '.pedsnp'
        copyfile(self.dataset.bim, pedsnp_filepath)
        return pedsnp_filepath

    def _create_pedind(self):
        ped = pd.read_table(self.dataset.ped, header=None, sep='\s+',
                            dtype=str)
        # ^ The 'dtype=str' is important to keep sample IDs and family IDs
        # untransformed. Otherwise, pandas tries to read them as numbers
        # and might change them.
        ped[5] = ped[5].replace(-9, 1)
        # ^ EIGENSOFT Q&A says we should replace 'funny values' like -9 for 1,
        # otherwise it sets 'ignore' to the samples with those values in the
        # 6th column.
        pedind = ped.ix[:, :5]  # .pedind = the first 6 columns of .ped
        pedind_filepath = self.dataset.label_path + '.pedind'
        pedind.to_csv(pedind_filepath, sep=' ', header=False, index=False)
        return pedind_filepath

    def _parse_evec_file(self, df):
        # The results file joins FID and IID with a ':'.
        df['FID'] = df[0].map(lambda s: s.split(':')[0])
        df['IID'] = df[0].map(lambda s: s.split(':')[1])
        df = df.set_index(['FID', 'IID'])
        df = df.drop(0, axis=1)
        df = df.ix[:, df.columns[:-1]]
        # ^ Removes a '???/Case/Control' col
        df.columns = ['PC{}'.format(column) for column in df.columns]
        df = df.join(self.dataset.samplegroup.samples).reset_index()
        df = df.set_index(['phenotype', 'sexcode', 'FID', 'IID'])
        df = df.sort_index()
        return df

    def _read_eval_file(self, eval_filename):
        df = pd.read_table(eval_filename, header=None, names=['variance'])
        df.index = ['PC{}'.format(ix + 1) for ix in df.index]
        df.index.name = 'component'
        df['ratio'] = df['variance'] / df['variance'].sum()
        df['percentage'] = df['ratio'].map(percentage_fmt)
        return df

    def _read_log(self):
        self._logfile
        return 'Not implemented yet :('
        # TODO read the logfile!
