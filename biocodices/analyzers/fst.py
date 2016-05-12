import pandas as pd
from helpers.plink import Plink


class Fst:
    @staticmethod
    def run(dataset, level):
        clusters_file = dataset.samplegroup.clusters_filepath(level)
        plink = Plink(dataset.bedfile)
        out_label = '{}.{}.{}'.format(dataset.samplegroup.label,
                                      dataset.panel.label, level)
        fst_file = plink.fst(clusters_file, out=out_label) + '.fst'
        fst_file_fields = ['chr', 'rs_id', 'position', 'NMISS', 'Fst']
        df = pd.read_table(fst_file, index_col='rs_id', skiprows=1,
                           names=fst_file_fields)
        return df
