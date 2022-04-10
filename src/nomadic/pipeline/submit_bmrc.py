import os
import re
import uuid
import configparser
from dataclasses import dataclass

from nomadic.lib.parsing import build_parameter_dict
from nomadic.lib.generic import print_header, print_footer


SUBMISSION_ID = str(uuid.uuid4())[:5]
BMRC_TEMPLATE = "bmrcs/template/bmrc_universal-template.sh"
SCRIPT_PATH = "configs/bmrc_scripts.ini"
PIPELINE_PATH = "configs/bmrc_pipeline.ini"


# ================================================================================
# Create BMRC shell scripts and qsub statements for a given script
#
# ================================================================================


class BmrcScriptBuilder:

    submission_id = SUBMISSION_ID
    bmrc_template = BMRC_TEMPLATE

    def __init__(self, script, job_name, queue="short.qc"):
        """
        Class to create BMRC shell scripts and qsub statements for
        a given input `script`

        params
            script: str
                How the script is called from the command line,
                e.g. `python run_mapping.py` or `nomadic map`, &c.
            job_name: str
                Job name for the script on BMRC computer cluster.
            queue: str
                Queue to submit job to on BMRC. By default sent
                to short queue.

        """

        self.script = script
        self.job_name = f"{job_name}-{self.submission_id}"
        self.log_name = job_name
        self.queue = queue
        self._create_log_directories()
        self.bmrc_output_path = None

    def _create_log_directories(self, log_dir="logs"):
        """Create the log directories if they do not exist"""
        streams = ["stdout", "stderr"]
        for stream in streams:
            stream_dir = f"{log_dir}/{stream}_{self.log_name}"
            if not os.path.isdir(stream_dir):
                print(f"Creating log directory: {stream_dir}")
                os.makedirs(stream_dir)

    def create_bmrc_script(
        self, kwargs, start_barcode=None, end_barcode=None, output_dir="bmrcs"
    ):
        """Generate a shell script for submission to BMRC"""

        extra_qsub_args = ""
        extra_script_args = ""

        # Check if it is a GPU submission
        if self.queue.endswith(".qg"):
            extra_qsub_args += "#$ -l gpu=1\n"

        # Prepare if array job
        if start_barcode is not None and end_barcode is not None:
            extra_qsub_args += f"#$ -t {start_barcode}-{end_barcode}:1"
            extra_script_args = "-b ${SGE_TASK_ID}"

        # Load template
        with open(self.bmrc_template, "r") as bmrc_template_file:

            # Replace variables in template
            bmrc_template_str = bmrc_template_file.read()
            output_str = bmrc_template_str.format(
                name=self.job_name,
                queue=self.queue,
                log=self.log_name,
                extra_qsub_args=extra_qsub_args,
                script=self.script,
                script_args=kwargs,  # " ".join([f"{k} {v}" for k, v in kwargs.items()]),
                extra_script_args=extra_script_args,
            )

            # Write output
            self.bmrc_output_path = f"{output_dir}/bmrc_{self.job_name}.sh"
            with open(self.bmrc_output_path, "w") as output:
                output.write(output_str)
            os.chmod(self.bmrc_output_path, 0o777)

        print(f"Submission script written to: {self.bmrc_output_path}")

    def get_qsub_statement(self, hold_for=None, array_hold=False):
        """
        Return a `qsub` statement for submission to BMRC

        """

        # Construct optional string indicating job dependency
        if hold_for:
            hold_str = " -hold_jid"
            hold_str += "_ad" if array_hold else ""
            hold_str += f" {hold_for}"
        else:
            hold_str = ""

        statement = f"qsub{hold_str} {self.bmrc_output_path}\n"

        return statement


# ================================================================================
# Define BMRC Script Factories for various scripts in NOMADIC pipeline
#
# ================================================================================


def create_bmrc_script_dict(script_path):
    """Create a dictionary of BMRC scripts"""

    # Create ConfigParser objects
    config = configparser.ConfigParser()
    config.read(script_path)

    # Iterate over scripts, create BMRC script objects
    scripts = {}
    for section in config.sections():
        scripts[section] = BmrcScriptBuilder(
            script=config.get(section, "script"),
            job_name=config.get(section, "job_name"),
            queue=config.get(section, "queue", fallback="short.qc"),
        )

    return scripts


# ================================================================================
# Define the BMRC Pipeline as a collection of (optionally) dependent scripts
#
# ================================================================================


@dataclass
class Script:
    builder: BmrcScriptBuilder
    args: dict
    array_job: bool = False
    dependency: str = None
    array_dependency: bool = False


def format_arguments(argument_template, params):
    """Format an argument string based on the parameter dictionary"""

    # Define an argument mapping
    # - Probably should either not exist, or exist elsewhere
    argument_map = {
        "-e": params["expt_dir"],
        "-c": params["config"],
        "-m": "hac",
        "-k": "native",
    }

    flags = re.findall("-[a-z]", argument_template)
    values = [argument_map[flag] for flag in flags]

    return argument_template.format(*values)


def create_bmrc_pipeline(pipeline_path, params, scripts):
    """Create a BMRC submission pipeline"""

    # Prepare ConfigParser instance
    config = configparser.ConfigParser()
    config.read(pipeline_path)

    # Create pipeline list
    pipeline = []
    for section in config.sections():
        pipeline_script = Script(
            scripts[section],
            args=format_arguments(config.get(section, "arguments"), params),
            array_job=config.getboolean(section, "array_job", fallback=None),
            dependency=scripts[config.get(section, "dependency")]
            if config.has_option(section, "dependency")
            else None,
            array_dependency=config.getboolean(
                section, "array_dependency", fallback=None
            ),
        )
        pipeline.append(pipeline_script)

    return pipeline


# ================================================================================
# Main script, run from `cli.py`
#
# ================================================================================


def main(expt_dir, config, barcode):  # barcode we don't need
    # PARSE INPUTS
    script_descrip = "NOMADIC: Prepare BMRC Pipeline Submission"
    t0 = print_header(script_descrip)
    params = build_parameter_dict(expt_dir, config, barcode)

    # PREPARE PIPELINE
    scripts = create_bmrc_script_dict(SCRIPT_PATH)
    pipeline = create_bmrc_pipeline(
        pipeline_path=PIPELINE_PATH, params=params, scripts=scripts
    )

    # Get start and end barcode, for array jobs
    start_barcode = int(min(params["barcodes"])[-2:])
    end_barcode = int(max(params["barcodes"])[-2:])

    # CREATE SUBMISSION FILE
    submission_fn = f"submit_pipeline-{SUBMISSION_ID}.sh"
    with open(submission_fn, "w") as fn:

        for script in pipeline:

            # Create the BMRC script
            script.builder.create_bmrc_script(
                kwargs=script.args,
                start_barcode=start_barcode if script.array_job else None,
                end_barcode=end_barcode if script.array_job else None,
            )

            # Write qsub
            qsub_statement = script.builder.get_qsub_statement(
                hold_for=script.dependency.job_name if script.dependency else None,
                array_hold=script.array_dependency,
            )
            fn.write(qsub_statement)
    os.chmod(submission_fn, 0o777)
    print(f"Submission file written to: {submission_fn}")
    print_footer(t0)
