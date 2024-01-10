#!/usr/bin/env python3

import argparse
import pathlib
import os
import re
import pandas as pd
import subprocess
import hashlib


def remove_dna_symbols(string):
    pattern = r"[ACGT]{5,}"
    result = re.sub(pattern, "", string)
    return result


def remove_lane_suffix(string):
    pattern = r"_L0*\d*"
    result = re.sub(pattern, "", string)
    return result


def sort_list_of_files(files_list):
    def natural_key(string):
        return [int(s) if s.isdigit() else s for s in re.split(r"(\d+)", string)]

    return sorted(files_list, key=lambda x: natural_key(x))


def validate_fastq_columns(fastq_1, fastq_2):
    pattern = re.compile(r"(.*)([R]?[12])\.([fastq]+)\.gz$")

    for fq1, fq2 in zip(fastq_1, fastq_2):
        match1 = pattern.match(fq1)
        match2 = pattern.match(fq2)

        if match1 and match2:
            base1 = match1.group(1)
            base2 = match2.group(1)
            if base1 != base2:
                raise ValueError(f"Failed to match fq1 and fq2 files:\n{fq1}\n{fq2}")
        else:
            raise ValueError(f"Failed to match fq1 and fq2 files:\n{fq1}\n{fq2}")


def validate_samplesheet(pd_df):
    if len(pd.unique(pd_df["sample"])) != len(pd_df["sample"]):
        raise ValueError(f"Failed to create a sample-sheet with unique sample IDs")

    validate_fastq_columns(list(pd_df["fastq_1"]), list(pd_df["fastq_2"]))


def compute_md5(file_path):
    with open(file_path, "rb") as file:
        data = file.read(10_000_000)
        md5_hash = hashlib.md5(data).hexdigest()

    print(f"Sample: {os.path.basename(file_path)}")

    return md5_hash


def remove_duplicate_samples(pd_df):
    print("Checking for duplicate samples and removing them...")

    pd_df["md5sum"] = pd_df.apply(
        lambda row: compute_md5(row["fastq_1"]) + compute_md5(row["fastq_2"]), axis=1
    )

    duplicates = pd_df.duplicated(subset="md5sum")

    if duplicates.any():
        print("Following samples have duplicates:")
        print(pd_df[duplicates].iloc[:, 0])
    else:
        print("No duplicates detected...")

    df_unique = pd_df.drop_duplicates(subset="md5sum", keep="first")
    return df_unique.iloc[:, :4]


def extract_r1_r2_files(list_of_files):
    list_of_files_sorted = sort_list_of_files(list_of_files)
    list_of_files_R1 = [
        f
        for f in list_of_files_sorted
        if len(re.findall(r"_[R]?1\.([fastq]+)\.gz", str(f))) > 0
    ]
    list_of_files_R2 = [
        f
        for f in list_of_files_sorted
        if len(re.findall(r"_[R]?2\.([fastq]+)\.gz", str(f))) > 0
    ]

    if len(list_of_files_R1) != len(list_of_files_R2):
        raise ValueError("Number of R1 and R2 files do not match")

    return list_of_files_R1, list_of_files_R2


def get_common_literals(list_of_lists):
    if len(list_of_lists) == 0:
        return []

    common_elements = set(list_of_lists[0])

    for sublist in list_of_lists[1:]:
        common_elements = common_elements.intersection(sublist)

    common_elements = list(common_elements)

    return sorted(common_elements)


def get_unique_elements(input_list):
    unique_elements = []
    seen_elements = set()
    for item in input_list:
        if item not in seen_elements:
            unique_elements.append(item)
            seen_elements.add(item)
    return unique_elements


def create_sample_ids_from_files_list(list_of_files_R1):
    file_names_normalized = [
        f.replace(".fastq.gz", "")
        .replace(".fq.gz", "")
        .replace("-", "_")
        .replace("/", "_")
        .replace("R1", "")
        for f in list_of_files_R1
    ]

    file_name_literals = [
        [l for l in f.split("_") if l != ""] for f in file_names_normalized
    ]

    common_literals = get_common_literals(file_name_literals)

    cleaved_names = []
    for f in file_names_normalized:
        cleaved_name = f
        for l in common_literals:
            cleaved_name = cleaved_name.replace(f"_{l}_", "__")

        cleaved_name = remove_dna_symbols(cleaved_name)

        cleaved_name_literals = [
            e for e in get_unique_elements(cleaved_name.split("_")) if e != ""
        ]

        cleaved_name = ""
        is_first = True
        for em in cleaved_name_literals:
            cleaved_name += em if is_first else "_" + em
            is_first = False

        cleaved_names.append(cleaved_name)

    sample_ids = cleaved_names

    return sample_ids


def save_samplesheet(exp_name, list_of_files_R1, list_of_files_R2, sample_ids):
    strandedness = ["auto" for _ in sample_ids]
    file_data = pd.DataFrame(
        {
            "sample": sample_ids,
            "fastq_1": list_of_files_R1,
            "fastq_2": list_of_files_R2,
            "strandedness": strandedness,
        }
    )

    validate_samplesheet(file_data)
    file_data.sort_values(by=["sample"], inplace=True)

    file_data_dedup = remove_duplicate_samples(file_data)

    file_data_dedup["sample"] = file_data_dedup["sample"].apply(remove_lane_suffix)

    file_data_dedup.to_csv(f"{exp_name}_samplesheet.csv", index=False)


def make_samplesheet_from_metadata_file(file_path, exp_name):
    file_data = pd.read_excel(file_path, sheet_name="Samplesheet")

    sample_id_col = "isolate"

    file_data.loc[:, "sample"] = file_data[sample_id_col]
    file_data.loc[:, "fastq_1"] = (
        file_data["directory"] + "/" + file_data["file_name_F"]
    )
    file_data.loc[:, "fastq_2"] = (
        file_data["directory"] + "/" + file_data["file_name_R"]
    )
    file_data.loc[:, "strandedness"] = "auto"

    file_data = file_data[["sample", "fastq_1", "fastq_2", "strandedness"]]

    validate_samplesheet(file_data)
    file_data.sort_values(by=["sample"], inplace=True)

    file_data_dedup = remove_duplicate_samples(file_data)
    file_data_dedup.to_csv(f"{exp_name}_samplesheet.csv", index=False)


def make_samplesheet_from_folder(file_path, exp_name):
    if os.path.isfile(file_path):
        raise ValueError(
            "The provided path is for a file. Path to an input folder is required"
        )

    fastq_gz_list = [str(f) for f in file_path.glob("*.fastq.gz")]
    fq_gz_list = [str(f) for f in file_path.glob("*.fq.gz")]

    list_of_files = fastq_gz_list + fq_gz_list

    if len(list_of_files) < 1:
        raise ValueError(
            "Could not find any fastq.gz or fq.gz files in the command output"
        )

    list_of_files_R1, list_of_files_R2 = extract_r1_r2_files(list_of_files)
    sample_ids = create_sample_ids_from_files_list(list_of_files_R1)
    save_samplesheet(exp_name, list_of_files_R1, list_of_files_R2, sample_ids)


def make_samplesheet_from_command(input_path_or_command, exp_name):
    result = subprocess.run(
        input_path_or_command, shell=True, capture_output=True, text=True
    )

    if result.returncode != 0:
        raise ValueError(f"Failed to execute the provided command...\n{result.stderr}")

    list_of_files = [
        f
        for f in result.stdout.split("\n")
        if f != "" and (f.endswith(".fq.gz") or f.endswith(".fastq.gz"))
    ]

    if len(list_of_files) < 1:
        raise ValueError(
            "Could not find any fastq.gz or fq.gz files in the command output"
        )

    list_of_files_R1, list_of_files_R2 = extract_r1_r2_files(list_of_files)
    sample_ids = create_sample_ids_from_files_list(list_of_files_R1)
    save_samplesheet(exp_name, list_of_files_R1, list_of_files_R2, sample_ids)

def main():
    parser = argparse.ArgumentParser(
        prog="make-sample-sheet",
        description="Read an RNASeq input folder or metadata file and create a sample-sheet.csv",
    )
    parser.add_argument(
        "path",
        help="RNASeq input folder path or metadata file path or a bash command which lists all the fastq samples",
    )
    parser.add_argument(
        "experiment-name",
        help="RNASeq experiment name",
    )
    parser.add_argument("-v", action="version", version="%(prog)s v0.3")

    args = vars(parser.parse_args())

    input_path_or_command = args["path"]
    exp_name = args["experiment-name"]

    print("Creating sample sheet...")

    if os.path.isfile(input_path_or_command):
        make_samplesheet_from_metadata_file(input_path_or_command, exp_name)
    elif os.path.isdir(input_path_or_command):
        make_samplesheet_from_folder(pathlib.Path(input_path_or_command), exp_name)
    else:
        make_samplesheet_from_command(input_path_or_command, exp_name)

if __name__ == "__main__":
    main()
