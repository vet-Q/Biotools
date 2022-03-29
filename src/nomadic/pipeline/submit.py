import os
import sys
import uuid
from dataclasses import dataclass
from nomadic.lib.parsing import basic_input_parser


SUBMISSION_ID = str(uuid.uuid4())[:5]
BMRC_TEMPLATE = "bmrcs/template/bmrc_universal-template.sh"


class BmrcScriptBuilder:

    submission_id = SUBMISSION_ID
    bmrc_template = BMRC_TEMPLATE

    def __init__(self, script, job_name, queue="short.qc"):
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
        """Generate a BMRC submission script"""

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
                script_args=" ".join([f"{k} {v}" for k, v in kwargs.items()]),
                extra_script_args=extra_script_args,
            )

            # Write output
            self.bmrc_output_path = f"{output_dir}/bmrc_{self.job_name}.sh"
            with open(self.bmrc_output_path, "w") as output:
                output.write(output_str)
            os.chmod(self.bmrc_output_path, 0o777)

        print(f"Submission script written to: {self.bmrc_output_path}")

    def get_qsub_statement(self, hold_for=None, array_hold=False):
        """Get a qsub statement for the script"""

        # Construct optional string indicating job dependency
        if hold_for:
            hold_str = " -hold_jid"
            hold_str += "_ad" if array_hold else ""
            hold_str += f" {hold_for}"
        else:
            hold_str = ""

        statement = f"qsub{hold_str} {self.bmrc_output_path}\n"

        return statement


@dataclass
class Script:
    builder: BmrcScriptBuilder
    args: dict
    array_job: bool = False
    dependency: str = None
    array_dependency: bool = False


# Define types of BMRC Scripts
Basecalling = BmrcScriptBuilder(
    script="./scripts/run_guppy-basecall.sh", job_name="basecall", queue="short.qg"
)
Barcoding = BmrcScriptBuilder(
    script="./scripts/run_guppy-barcode.sh", job_name="barcode", queue="short.qg"
)
Mapping = BmrcScriptBuilder(script="python run_mapping-update.py", job_name="map")
Remapping = BmrcScriptBuilder(script="python run_mapping-hs.py", job_name="remap")
QCBams = BmrcScriptBuilder(script="python run_qc-bam-update.py", job_name="qcbam")


def main():
    # PARSE INPUTS
    script_descrip = "NOMADIC: Prepare BMRC Pipeline Submission"
    params = basic_input_parser(sys.argv, descrip=script_descrip)

    # Get start and end barcode, for array jobs
    start_barcode = int(min(params["barcodes"])[-2:])
    end_barcode = int(max(params["barcodes"])[-2:])

    # Define arguments
    basic_args = {
        "-e": params["expt_dir"],
        "-c": params["config"],
    }
    basecall_args = {
        "-e": params["expt_dir"],
        "-m": "hac",
    }
    barcode_args = {
        "-e": params["expt_dir"],
        "-m": "hac",
        "-k": "native",
    }

    # Scripts to run
    scripts = [
        Script(Basecalling, basecall_args),
        Script(Barcoding, barcode_args, dependency=Basecalling),
        Script(Mapping, basic_args, array_job=True, dependency=Barcoding),
        Script(
            Remapping,
            basic_args,
            array_job=True,
            dependency=Mapping,
            array_dependency=True,
        ),
        Script(
            QCBams,
            basic_args,
            array_job=True,
            dependency=Remapping,
            array_dependency=True,
        ),
    ]

    # CREATE SUBMISSION FILE
    submission_fn = f"submit_pipeline-{SUBMISSION_ID}.sh"
    with open(submission_fn, "w") as fn:

        for script in scripts:

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


if __name__ == "__main__":
    main()
