import subprocess
from datetime import datetime


class ProgramCaller:
    def __init__(self, command, log_filepath):
        self.command = command
        self.log_filepath = log_filepath

    def run(self):
        self.t1 = datetime.now()

        with open(self.log_filepath, 'w+') as log_file:
            self._write_command(log_file)

        with open(self.log_filepath, 'a+') as log_file:
            self.proc = subprocess.run(self.command.split(' '), stdout=log_file,
                                       stderr=subprocess.STDOUT)
            self.t2 = datetime.now()
            self._write_timestamp(log_file)
        return self.proc

    def _write_command(self, file_handle):
        tmpl = 'Started at {}:\n{}\n---\n\n'
        file_handle.write(tmpl.format(self._timestamp(), self.command))

    def _write_timestamp(self, file_handle):
        elapsed = (self.t2 - self.t1).seconds
        tmpl = '\n---\nFinished at {}.\nTook {} seconds.\n'
        file_handle.write(tmpl.format(self._timestamp(), elapsed))

    def _timestamp(self):
        return datetime.now().strftime('%Y-%b-%d %X')
