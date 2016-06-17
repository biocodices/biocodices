from os.path import join
from .config import Config


class Resource:
    @staticmethod
    def available_resources():
        return Config('resources')

    def __new__(self, label):
        """
        Handles resources filepaths. You can query a deep key from the yaml
        by separating nested keys with a ':', like:
            Resource('top_key:deep_key')
        """
        resources = Config('resources')
        base_dir = resources['base_dir']

        # Hack to do the nested lookups in the dict:
        value = resources  # Top level
        for key in label.split(':'):
            value = value[key]  # Gets one level deeper each time

        return join(base_dir, value)
