import yaml
from os.path import join, expanduser


class Config:
    def __new__(self, config_label):
        """
        Expects a ~/.biocodices dir with yml config files in it.
        Call directly Config('<config_label'>) to get a dictionary with settings
        """
        self.base_dir = expanduser('~/.biocodices')
        return yaml.load(open(join(self.base_dir, config_label + '.yml')))
