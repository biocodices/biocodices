import yaml

from os.path import join, expanduser


class Config:
    BASE_DIR = expanduser('~/biocodices/settings')

    def __new__(self, config_label):
        return yaml.load(open(join(self.BASE_DIR, config_label + '.yml')))
