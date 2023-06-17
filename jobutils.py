# Copyright (c) 2023 The Regents of the Huazhong University of Science and Technology
# All rights reserved.
# This software is licensed under the BSD 3-Clause License.
# Authors: Teng Zhang (2022010236@hust.edu.cn)
"""
This module adds jobcontrol related command line flags and job utilities.
"""
import json
from signac import job
from nemd import environutils
import os
import argparse


RUN_NEMD = 'run_nemd'
OUTFILE = 'outfile'
LOGFILE = 'logfile'
OUTFILES = 'outfiles'
TARGS = 'targs'
FLAG_INTERACTIVE = '-INTERACTIVE'
FLAG_JOBNAME = '-JOBNAME'
FLAG_DEBUG = '-DEBUG'
FLAG_SEED = '-seed'
FLAG_CLEAN = '-clean'
FLAG_JTYPE = '-jtype'
FLAG_CPU = '-cpu'
FLAG_PRJ_PATH = '-prj_path'
PREREQ = 'prereq'

FINISHED = 'Finished.'
FILE = "$FILE"
ARGS = 'args'
KNOWN_ARGS = 'known_args'
UNKNOWN_ARGS = 'unknown_args'
FN_DOCUMENT = job.Job.FN_DOCUMENT
TASK = 'task'
AGGREGATOR = 'aggregator'

def type_positive_int(arg):
    value = type_int(arg)
    if value < 1:
        raise argparse.ArgumentTypeError(f'{value} is not a positive integer.')
    return value

def add_job_arguments(parser, arg_flags=None, jobname=None):
    """
    Add job control related flags.

    :param parser: the parser to add arguments
    :type parser: 'argparse.ArgumentParser'
    :param arg_flags: specific job control related flags to add
    :type arg_flags: list
    :param jobname: the default jobname
    :type jobname: str
    """
    if arg_flags is None:
        arg_flags = [FLAG_INTERACTIVE, FLAG_JOBNAME, FLAG_DEBUG, FLAG_CPU]
    # Workflow drivers may add the job control options a few times
    arg_flags = [
        x for x in arg_flags if x not in parser._option_string_actions
    ]
    if FLAG_INTERACTIVE in arg_flags:
        parser.add_argument(FLAG_INTERACTIVE,
                            dest=FLAG_INTERACTIVE[1:].lower(),
                            action='store_true',
                            help='Enable interactive mode')
    if FLAG_JOBNAME in arg_flags:
        parser.add_argument(
            FLAG_JOBNAME,
            dest=FLAG_JOBNAME[1:].lower(),
            default=jobname,
            help='The jobnamee based on which filenames are created.')
    if FLAG_DEBUG in arg_flags:
        parser.add_argument(
            FLAG_DEBUG,
            action='store_true',
            dest=FLAG_DEBUG[1:].lower(),
            help='Enable debug mode (e.g. extra printing and files)')
    if FLAG_CPU in arg_flags:
        parser.add_argument(FLAG_CPU,
                            type=type_positive_int,
                            dest=FLAG_CPU[1:].lower(),
                            default=round(os.cpu_count() / 2),
                            help='Number of CPU processors.')


def get_arg(args, flag, val=None):
    """
    Get the value after the flag in command arg list.

    :param args: the arg list
    :type args: list
    :param flag: set the value after this flag
    :type flag: str
    :param val: the default value if the flag doesn't exist
    :type val: str
    :return: the value after the flag
    :rtype: str
    """
    try:
        idx = args.index(flag)
    except ValueError:
        return val
    else:
        return args[idx + 1]


def set_arg(args, flag, val):
    """
    Set the value after the flag in command arg list.

    :param args: the arg list
    :type args: list
    :param flag: set the value after this flag
    :type flag: str
    :param val: the new value
    :type val: str
    :return: the modified arg list
    :rtype: list
    """
    if flag not in args:
        args.extend([flag, ''])
    idx = args.index(flag)
    args[idx + 1] = val
    return args


def add_outfile(outfile,
                jobname=None,
                default_jobname=None,
                job=None,
                document=FN_DOCUMENT,
                set_file=False,
                log_file=False):
    """
    Register the outfile to the job control.

    :param outfile: the outfile to be registered
    :type outfile: str
    :param jobname: register the file under this jobname
    :type jobname: str
    :param default_jobname: use this jobname if jobname is undefined and cannot
        be found from the environment
    :type default_jobname: str
    :param job: register outfile to this job document
    :type job: 'signac.contrib.job.Job'
    :param document: the job control information is saved into this file if job
        not provided
    :type document: str
    :param set_file: set this file as the single output file
    :type set_file: bool
    :param log_file: set this file as the log file
    :type log_file: bool
    """
    if jobname is None:
        jobname = environutils.get_jobname(default_jobname)
    if job:
        data = job.document
    else:
        try:
            with open(document) as fh:
                data = json.load(fh)
        except FileNotFoundError:
            data = {}
    data.setdefault(OUTFILES, {})
    data[OUTFILES].setdefault(jobname, [])
    if outfile not in data[OUTFILES][jobname]:
        data[OUTFILES][jobname].append(outfile)
    if set_file:
        data.setdefault(OUTFILE, {})
        data[OUTFILE].setdefault(jobname, outfile)
    if log_file:
        data.setdefault(LOGFILE, {})
        data[LOGFILE].setdefault(jobname, outfile)
    if job:
        return
    with open(document, 'w') as fh:
        json.dump(data, fh)
