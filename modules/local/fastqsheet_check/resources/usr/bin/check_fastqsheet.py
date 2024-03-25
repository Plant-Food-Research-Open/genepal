#!/usr/bin/env python3

import os
import sys
import errno
import argparse

# https://github.com/nf-core/rnaseq
# MIT: https://github.com/nf-core/rnaseq/blob/master/LICENSE
#
# Changes:
#
# 1. Formatted with black
# 2. Added checks for the fifth column: target_assemblies
# 3. Removed strandedness


def parse_args(args=None):
    Description = (
        "Reformat nf-core/rnaseq style samplesheet file and check its contents."
    )
    Epilog = 'Example usage: python check_samplesheet.py <FILE_IN> "target_assembly_a,target_assembly_b" <FILE_OUT>'

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input samplesheet file.")
    parser.add_argument("TARGET_ASSEMBLIES", help="Permissible target assemblies")
    parser.add_argument("FILE_OUT", help="Output file.")
    return parser.parse_args(args)


def make_dir(path):
    if len(path) > 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise exception


def print_error(error, context="Line", context_str=""):
    error_str = f"ERROR: Please check samplesheet -> {error}"
    if context != "" and context_str != "":
        error_str = f"ERROR: Please check samplesheet -> {error}\n{context.strip()}: '{context_str.strip()}'"
    print(error_str)
    sys.exit(1)


def check_samplesheet(file_in, file_out, permissible_target_assemblies):
    """
    This function checks that the samplesheet follows the following structure:

    sample,fastq_1,fastq_2,target_assemblies
    SAMPLE_PE,SAMPLE_PE_RUN1_1.fastq.gz,SAMPLE_PE_RUN1_2.fastq.gz,red5;red3
    SAMPLE_PE,SAMPLE_PE_RUN2_1.fastq.gz,SAMPLE_PE_RUN2_2.fastq.gz,red5;red3
    SAMPLE_SE,SAMPLE_SE_RUN1_1.fastq.gz,,red5
    """

    sample_mapping_dict = {}
    with open(file_in, "r", encoding="utf-8-sig") as fin:
        ## Check header
        MIN_COLS = 4
        HEADER = ["sample", "fastq_1", "fastq_2", "target_assemblies"]
        header = [x.strip('"') for x in fin.readline().strip().split(",")]
        if header[: len(HEADER)] != HEADER:
            print(
                f"ERROR: Please check samplesheet header -> {','.join(header)} != {','.join(HEADER)}"
            )
            sys.exit(1)

        ## Check sample entries
        for line in fin:
            if line.strip():
                lspl = [x.strip().strip('"') for x in line.strip().split(",")]

                ## Check valid number of columns per row
                if len(lspl) < len(HEADER):
                    print_error(
                        f"Invalid number of columns (minimum = {len(HEADER)})!",
                        "Line",
                        line,
                    )

                num_cols = len([x for x in lspl[: len(HEADER)] if x])
                if num_cols < MIN_COLS:
                    print_error(
                        f"Invalid number of populated columns (minimum = {MIN_COLS})!",
                        "Line",
                        line,
                    )

                ## Check sample name entries
                sample, fastq_1, fastq_2, target_assemblies = lspl[
                    : len(HEADER)
                ]

                if sample.find(" ") != -1:
                    print(
                        f"WARNING: Spaces have been replaced by underscores for sample: {sample}"
                    )
                    sample = sample.replace(" ", "_")
                if not sample:
                    print_error("Sample entry has not been specified!", "Line", line)

                ## Check FastQ file extension
                for fastq in [fastq_1, fastq_2]:
                    if fastq:
                        if fastq.find(" ") != -1:
                            print_error("FastQ file contains spaces!", "Line", line)
                        if not fastq.endswith(".fastq.gz") and not fastq.endswith(
                            ".fq.gz"
                        ):
                            print_error(
                                "FastQ file does not have extension '.fastq.gz' or '.fq.gz'!",
                                "Line",
                                line,
                            )

                ## Auto-detect paired-end/single-end
                sample_info = []  ## [single_end, fastq_1, fastq_2]
                if sample and fastq_1 and fastq_2:  ## Paired-end short reads
                    sample_info = ["0", fastq_1, fastq_2]
                elif sample and fastq_1 and not fastq_2:  ## Single-end short reads
                    sample_info = ["1", fastq_1, fastq_2]
                else:
                    print_error(
                        "Invalid combination of columns provided!", "Line", line
                    )

                ## Check if the target assemblies are permissible
                target_assemblies_list = sorted(
                    [x.strip() for x in target_assemblies.strip().split(";")]
                )

                for assembly in target_assemblies_list:
                    if assembly in permissible_target_assemblies:
                        continue

                    print_error(
                        f"Target assembly '{assembly}' is not one of {permissible_target_assemblies}!",
                        "Line",
                        line,
                    )

                ## Create sample mapping dictionary = {sample: [[ single_end, fastq_1, fastq_2, strandedness, [target_assemblies] ]]}
                sample_target_assemblies = ";".join(target_assemblies_list)
                sample_info = (
                    sample_info + lspl[len(HEADER) :] + [sample_target_assemblies]
                )
                if sample not in sample_mapping_dict:
                    sample_mapping_dict[sample] = [sample_info]
                else:
                    if sample_info in sample_mapping_dict[sample]:
                        print_error(
                            "Samplesheet contains duplicate rows!", "Line", line
                        )
                    else:
                        sample_mapping_dict[sample].append(sample_info)

    ## Write validated samplesheet with appropriate columns
    if len(sample_mapping_dict) > 0:
        out_dir = os.path.dirname(file_out)
        make_dir(out_dir)
        with open(file_out, "w") as fout:
            fout.write(
                ",".join(
                    [
                        "sample",
                        "single_end",
                        "fastq_1",
                        "fastq_2",
                        "target_assemblies",
                    ]
                    + header[len(HEADER) :]
                )
                + "\n"
            )
            for sample in sorted(sample_mapping_dict.keys()):
                ## Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
                if not all(
                    x[0] == sample_mapping_dict[sample][0][0]
                    for x in sample_mapping_dict[sample]
                ):
                    print_error(
                        f"Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end!",
                        "Sample",
                        sample,
                    )

                ## Check that multiple runs of the same sample have same target assemblies
                if not all(
                    x[3] == sample_mapping_dict[sample][0][3]
                    for x in sample_mapping_dict[sample]
                ):
                    print_error(
                        f"Multiple runs of a sample must have the same target assemblies!",
                        "Sample",
                        sample,
                    )

                for idx, val in enumerate(sample_mapping_dict[sample]):
                    fout.write(",".join([f"{sample}_T{idx+1}"] + val) + "\n")
    else:
        print_error(f"No entries to process!", "Samplesheet: {file_in}")


def main(args=None):
    args = parse_args(args)
    permissible_target_assemblies = [
        x.strip() for x in args.TARGET_ASSEMBLIES.strip().split(",")
    ]
    check_samplesheet(args.FILE_IN, args.FILE_OUT, permissible_target_assemblies)


if __name__ == "__main__":
    sys.exit(main())
