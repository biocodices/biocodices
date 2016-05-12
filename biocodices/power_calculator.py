#!/usr/bin/env python3

import requests
import pandas as pd


class PowerCalculator:
    """
    Visits 'http://pngu.mgh.harvard.edu/~purcell/cgi-bin/cc2k.cgi' to get the
    number of cases needed for a case-control test. Checkout the parameters.
    """
    def case_control_for_discrete_trait(self, high_risk_allele_freq,
            prevalence, rel_risk_Aa, rel_risk_AA, d_prime, marker_allele_freq,
            number_of_cases, case_control_ratio=1, unselected_controls=True,
            typeI_error_rate=0.05, power=0.80):

        url = 'http://pngu.mgh.harvard.edu/~purcell/cgi-bin/cc2k.cgi'
        payload = {
            'fA': high_risk_allele_freq,
            'k': prevalence,
            'rAa': rel_risk_Aa,
            'rAA': rel_risk_AA,
            'dprime': d_prime,
            'm1': marker_allele_freq,
            'n': number_of_cases,
            'ccratio': case_control_ratio,
            'alpha': typeI_error_rate,
            'power': power,
        }

        response = requests.post(url, data=payload)
        tables = pd.read_html(response.text, header=0)
        # Last table on the page is the interesing one. No "id" on the <table>.
        allelic_1df_test = tables[-1]
        allelic_1df_test.columns = ['Alpha', 'Power', 'Ncases']
        allelic_1df_test.set_index('Alpha', inplace=True)
        allelic_1df_test.drop_duplicates(inplace=True)
        allelic_1df_test['Ncases'] = allelic_1df_test['Ncases'].astype('int32')

        return allelic_1df_test.loc[typeI_error_rate, 'Ncases']
