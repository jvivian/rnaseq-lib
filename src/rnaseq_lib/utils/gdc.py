#!/usr/bin/python
# -*- coding:utf-8 -*-
################################################
# File Name: mapFileUUID2submitterID.py
# Author: C.J. Liu
# Edited: John Vivian (8-30-2017)
# Mail: samliu@hust.edu.cn
# Created Time: Tue 15 Nov 2016 03:53:10 PM CST
################################################

import argparse
import os

import requests


def usage():
    description = '''
    Problem: 
        I choose project BRCA on GDC(https://gdc-portal.nci.nih.gov/projects/t) 
        and forward to WXS(1,050 cases, 2,175 bam files). Download manifest file from the summary tab. 
        In fact the manifest only contains the file_id, file_name and md5. 
        I can't tell sample information like TCGA barcode. 
        On the Files tab in the search result of GDC-Portal, 
        we can download sample information in form of json in some condition. 
        However, I can't get sample information from the BRCA-WXS json file. Use GDC-API files endpoint.
    Task: 
        Map file_id to sample sample id through GDC-API
         (https://gdc-docs.nci.nih.gov/API/Users_Guide/Search_and_Retrieval/#filters-specifying-the-quer
    '''
    use = """%(prog)s -i <manifest.txt>"""
    parser = argparse.ArgumentParser(description=description, usage=use)

    parser.add_argument("-i", "--input", dest='manifest', type=str,
                        help="Input manifest file downloaded from gdc-portal or file contains first column as file_id",
                        required=True)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')
    parser.add_argument('-o', '--output', type=str, default='./', help='Output location for metadata.')

    args = parser.parse_args()
    return args


def load_manifest(manifest):
    file_ids = list()
    with open(manifest, 'r') as foo:
        for line in foo:
            arr = line.rstrip().split("\t")
            if len(arr[0]) != 36:
                continue
            file_ids.append(arr[0])
    return file_ids


def make_params(file_ids):
    params = {
        "filters": {
            "op": "in",
            "content": {
                "field": "files.file_id",
                "value": file_ids
            }
        },
        "format": "TSV",
        "fields": "analysis.metadata.read_groups.is_paired_end,cases.samples.portions.analytes.aliquots.submitter_id,"
                  "cases.samples.sample_type",
        # The must be no space after comma
        "size": len(file_ids)
    }
    return params


def gdc_api(file_ids, manifest):
    output = manifest + '.map2submitterID'
    outputh = open(output, 'w')
    files_endpt = "https://gdc-api.nci.nih.gov/files"
    params = make_params(file_ids)
    response = requests.post(files_endpt, json=params)
    outputh.write(response.text)
    outputh.close()
    # print(response.text)
    return response.text


def metadata_from_manifest(manifest, output):
    file_ids = load_manifest(manifest)
    response = gdc_api(file_ids, manifest)
    with open(output, 'w') as f:
        f.write(response)


def main():
    args = usage()
    assert os.path.isdir(args.output)
    output = os.path.join(args.output, 'metadata.tsv')
    assert not os.path.exists(output)
    metadata_from_manifest(args.manifest, output)


if __name__ == '__main__':
    main()
