"""Create a dataset from raw data."""

from pathlib import Path
import re
from typing import Tuple
from xml.etree import ElementTree as ET

import numpy as np
import pandas as pd

KEEP_FIELDNAMES = (
    "disease_code",
    "tumor_tissue_site",
    "days_to_birth",
    "days_to_death",
    "days_to_last_followup",
    "vital_status",
    "clinical_T",
    "pathologic_T",
    "clinical_stage",
    "pathologic_stage",
    "measure_of_response",
    "primary_therapy_outcome_success",
    "days_to_drug_therapy_start",
    "days_to_drug_therapy_end",
    "person_neoplasm_cancer_status",
    "drug_name",
    "neoplasm_histologic_grade",
)


def clean_clinical_data(
    data: dict, keep_field_names: Tuple[str, ...] = KEEP_FIELDNAMES
):
    """Clean clinical data by removing empty fields as well as fields that won't be used.

    Args:
        data (dict): Clinical data.
        keep_field_names (Tuple[str]): Field names to keep.
    """
    for k, v in list(data.items()):
        if isinstance(v, dict):
            clean_clinical_data(v, keep_field_names)
            if len(v) == 0:
                del data[k]
        elif isinstance(v, list):
            if all(isinstance(val, str) for val in v) and k not in keep_field_names:
                del data[k]
                continue
            indices_to_remove = []
            for i, v_i in enumerate(v):
                if isinstance(v_i, dict):
                    clean_clinical_data(v_i, keep_field_names)
                    if len(v_i) == 0:
                        indices_to_remove.append(i)
            for i in indices_to_remove[::-1]:
                del v[i]
            if len(v) == 0:
                del data[k]
        elif v == "" or v == {} or k not in keep_field_names:
            del data[k]


def flatten_dict(d: dict) -> dict:
    """Flatten a nested dictionary to a single level.
    Only the lowest level key/value pairs are kept.

    Example:
    >>> flatten_dict({'a': {'b': 1, 'c': 2}, 'd': 3})
    {'b': 1, 'c': 2, 'd': 3}
    """

    def get_leaves(d: dict):
        for k, v in d.items():
            if isinstance(v, dict):
                for k2, v2 in get_leaves(v):
                    yield (k2, v2)
            else:
                yield (k, v)

    return dict(get_leaves(d))


def number_nested_lists(data):
    if isinstance(data, dict):
        for k, v in data.items():
            if isinstance(v, (dict, list)):
                data[k] = number_nested_lists(v)
    elif isinstance(data, list):
        data = {i: number_nested_lists(v) for (i, v) in enumerate(data)}
    return data


def restructure_clinical_dict(clinical_dict: dict) -> dict:
    """Restructure the clinical dictionary to be more easily accessible
    and to have a consistent structure."""
    new_clinical_dict = {}
    for patient, patient_data in clinical_dict.items():
        if not ("admin" in patient_data and "patient" in patient_data):
            continue
        new_data = {}
        new_data.update(patient_data["admin"])
        old_data: dict = patient_data["patient"]
        followups = old_data.get("follow_ups", {}).get("follow_up", [])
        if isinstance(followups, dict):
            followups = [followups]
        drugs = old_data.get("drugs", {}).get("drug", [])
        if isinstance(drugs, dict):
            drugs = [drugs]
        new_data["follow_ups"] = followups
        new_data["drugs"] = drugs
        for i, treatment in enumerate(new_data["drugs"]):
            new_data["drugs"][i] = flatten_dict(treatment)
        for k, v in old_data.items():
            if k in ("follow_ups", "drugs"):
                continue
            if isinstance(v, dict):
                # note, this does not maintain parent keys, only the lowest level field names
                new_data.update(flatten_dict(v))
            elif isinstance(v, list):
                for i, item in enumerate(v):
                    if isinstance(item, dict):
                        v[i] = flatten_dict(item)
                new_data[k] = v
            else:
                new_data[k] = v
        new_clinical_dict[patient] = new_data
    return number_nested_lists(new_clinical_dict)


def remove_http_text(xml_tag: str):
    """Remove http text inside curly braces from an XML tag."""
    return re.sub(r"\{http.*?\}", "", xml_tag)


def xml_to_dict_rec(root: ET.Element):
    """Recursively convert an XML element to a dictionary."""
    data = {}
    for child in root:
        child_data = xml_to_dict_rec(child)
        if child_data is not None:
            child_tag = remove_http_text(child.tag)
            # Handle multiple child elements with the same tag
            if child_tag in data:
                if not isinstance(data[child_tag], list):
                    data[child_tag] = [data[child_tag]]
                data[child_tag].append(child_data)
            else:
                data[child_tag] = child_data

    if data:
        return data

    return root.text.strip() if root.text else None


def xml_to_dict(xml_path: str):
    """Convert an XML file to a dictionary.

    Args:
        xml_path (str): Path to an XML file.

    Returns:
        dict: Dictionary representation of the XML file.
    """
    tree = ET.parse(xml_path)
    root = tree.getroot()
    return xml_to_dict_rec(root)


def parse_clinical_data(
    clinical_dir: str,
    keep_field_names: Tuple[str, ...] = KEEP_FIELDNAMES,
    TCGA_only: bool = True,
) -> dict:
    """Parse clinical data from XML files into a dictionary.

    Args:
        clinical_dir (str): Path to a directory containing clinical data XML files.
        keep_field_names (Tuple[str]): Field names to keep.

    Returns:
        dict: Dictionary representation of the clinical data.
    """
    clinical_data = {}
    clinical_dir = Path(clinical_dir)
    for xml_path in clinical_dir.glob("*.xml"):
        patient_id = str(xml_path).split(".")[-2]
        if TCGA_only and not patient_id.startswith("TCGA"):
            continue
        clinical_data[patient_id] = xml_to_dict(xml_path)

    clean_clinical_data(clinical_data, keep_field_names)
    clinical_data = restructure_clinical_dict(clinical_data)

    return clinical_data


def make_dataset():
    """Runs data processing scripts to turn raw data from (../raw) into
    cleaned data ready to be analyzed (saved in ../processed).

    Processing steps:
    1. RNA-seq data:
        1.1. Match RNA-seq profiles with clinical data.
        1.2. Where there are multiple profiles for a single patient,
            retain the profile with the highest mean expression.
        1.3. Remove genes that are not either protein-coding or lncRNA.
    2. Clinical data:
        1.1. Parse XML files into a dictionary.
        1.2. Remove fields that won't be used, enforce data types, and rename fields.
    """
