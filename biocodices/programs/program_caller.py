import subprocess
import re
from datetime import datetime
from biocodices.helpers.language import seconds_to_hms_string


class ProgramCaller:
    def __init__(self, command):
        self.command = re.sub(' +', ' ', command)
        # Trim multiple spaces, so they don't create a 'phantom' empty argument.

    def run(self, stdout_sink=None, stderr_sink=None, log_filepath=None):
        """
        Specify a log_filepath if you want log info with timestamps to be
        written there. Otherwise, you can redirect untoucheed stdout and
        stderr to filepaths of your choice.
        """
        start_time = datetime.now()

        self.log_filepath = log_filepath or 'biocodices_output.log'
        if not self.log_filepath.endswith('.log'):
            self.log_filepath += '.log'

        with open(self.log_filepath, 'w') as log_file:
            # This file write needs to be separated from the following one.
            # Otherwise, the written lines will appear at the bottom of the file.
            self._log_executed_command(log_file)

        with open(self.log_filepath, 'a') as log_file:
            stdout_sink_IO, stderr_sink_IO = self._define_sinks(
                stdout_sink, stderr_sink, log_file)

            try:
                self.process = subprocess.run(
                    self._command_to_arglist(self.command),
                    stdout=stdout_sink_IO,
                    stderr=stderr_sink_IO,
                    check=True)
            except subprocess.CalledProcessError:
                msg = (
                    '* This command failed:\n{}\n'
                    '* With stdout_sink:\n{}\n'
                    '* With stderr_sink:\n{}\n'
                ).format(
                    self.command,
                    (stdout_sink or self.log_filepath),
                    (stderr_sink or self.log_filepath)
                )
                print(msg)
                log_file.write('---\n' + msg)
                raise
            finally:
                # The log_file will be closed automatically because of the
                # 'with' context. For the other file handles, be cautious:
                if stdout_sink:
                    stdout_sink_IO.close()
                if stderr_sink:
                    stderr_sink_IO.close()

            self._log_elapsed_time(log_file, start_time, datetime.now())

    @staticmethod
    def _command_to_arglist(command):
        # GATK parameters include filter expressions as strings that must be
        # passed as a single argument, so I can't just do a command.split(' ').
        # I have to take them out first and restore them to the argument list
        # afterwards. For instance, '--filterExpression "QD < 2.0"' should be
        # passed as two arguments, not four.
        string_arguments = re.findall(r'(".*?")', command)
        dummy_keys = ['#{}#'.format(i) for i in range(len(string_arguments))]
        dummy_to_arg = dict(zip(dummy_keys, string_arguments))
        for key, arg in dummy_to_arg.items():
            # Replace the quoted expression with a dummy key
            command = command.replace(arg, key)

        arglist = command.split(' ')
        for i, arg in enumerate(arglist):
            if arg in dummy_to_arg:  # if the arg is actually a dummy key
                arglist[i] = dummy_to_arg[arg].replace('"', '')
                # Removing the quotations is necessary! GATK filter expressions
                # will come here surrounded in double quotes, but to pass the
                # expression correctly to GATK via subprocess, we need to
                # remove them (otherwise they are taken to be part of the
                # expression itself, they break the VCF and the filter is not
                # applied).

        return arglist

    @staticmethod
    def _define_sinks(stdout_sink, stderr_sink, log_file):
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

        return stdout_sink_IO, stderr_sink_IO

    def _log_executed_command(self, file_handle):
        start_msg = 'Started at {}\n{}\n---\n\n'
        file_handle.write(start_msg.format(self._timestamp(), self.command))

    @classmethod
    def _log_elapsed_time(cls, file_handle, t1, t2):
        seconds = (t2 - t1).seconds
        end_msg = '\n---\nFinished at {}\nTook '.format(cls._timestamp())
        end_msg += seconds_to_hms_string(seconds) + '.\n'
        file_handle.write(end_msg)

    @staticmethod
    def _timestamp():
        return datetime.now().strftime('%Y-%b-%d %X')
