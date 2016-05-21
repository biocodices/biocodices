import subprocess
from datetime import datetime
from biocodices.helpers.language import seconds_to_hms_string


class ProgramCaller:
    def __init__(self, command):
        self.command = command.replace('  ', ' ')
        # Trim double spaces, so they don't create a 'phantom' empty argument.

    def run(self, stdout_sink=None, stderr_sink=None, log_filepath=None):
        """
        Specify a log_filepath if you want log info with timestamps to be
        written there. Otherwise, you can redirect untoucheed stdout and
        stderr to filepaths of your choice.
        """
        self.t1 = datetime.now()

        self.log_filepath = log_filepath or 'biocodices_output.log'
        if not self.log_filepath.endswith('.log'):
            self.log_filepath += '.log'

        with open(self.log_filepath, 'w') as log_file:
            # This file write neeeds to be separated from the following one.
            # Otherwise, the written lines will appear at the bottom of the file.
            self._log_executed_command(log_file)

        with open(self.log_filepath, 'a') as log_file:
            # The following logic can be summarized as:
            # - If the user specified a sink file for some ouptut channel, left
            # it untouched and redirect it there.
            # - If the user didn't specify a sink file for some output channel,
            # put that in the log_file.
            # By output channels I mean: STDOUT, STDERR or both.
            if stdout_sink and stderr_sink:
                stdout_sink_IO = open(stdout_sink, 'w')
                stderr_sink_IO = open(stderr_sink, 'w')
            elif stdout_sink and not stderr_sink:
                stdout_sink_IO = open(stdout_sink, 'w')
                stderr_sink_IO = log_file
            elif not stdout_sink and stderr_sink:
                stdout_sink_IO = log_file
                stderr_sink_IO = open(stderr_sink, 'w')
            elif not stdout_sink and not stderr_sink:
                stdout_sink_IO = log_file
                stderr_sink_IO = subprocess.STDOUT  # Both outputs to the log

            try:
                self.proc = subprocess.run(self.command.split(' '),
                                           stdout=stdout_sink_IO,
                                           stderr=stderr_sink_IO,
                                           check=True)
            except subprocess.CalledProcessError as error:
                print('* This command failed:\n{}\n'.format(self.command))
                print('* With stdout_sink:\n{}\n'.format(stdout_sink or self.log_filepath))
                print('* With stderr_sink:\n{}\n'.format(stderr_sink or self.log_filepath))
                raise(error)
            finally:
                # The log_file will be closed automatically because of the
                # 'with' context. For the other file handles, be cautious:
                if stdout_sink:
                    stdout_sink_IO.close()
                if stderr_sink:
                    stderr_sink_IO.close()

            self.t2 = datetime.now()
            self._log_elapsed_time(log_file)
        return self.proc

    def _log_executed_command(self, file_handle):
        start_msg = 'Started at {}\n{}\n---\n\n'
        file_handle.write(start_msg.format(self._timestamp(), self.command))

    def _log_elapsed_time(self, file_handle):
        seconds = (self.t2 - self.t1).seconds
        end_msg = '\n---\nFinished at {}\nTook '.format(self._timestamp())
        end_msg += seconds_to_hms_string(seconds) + '.\n'
        file_handle.write(end_msg)

    def _timestamp(self):
        return datetime.now().strftime('%Y-%b-%d %X')
