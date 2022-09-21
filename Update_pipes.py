#!/usr/bin/env python3.9

import os
import paramiko
import argparse
import sys

class MySFTPClient(paramiko.SFTPClient):
    def put_dir(self, source, target):
        ''' Uploads the contents of the source directory to the target path. The
            target directory needs to exists. All subdirectories in source are 
            created under target.
        '''
        for item in os.listdir(source):
            if os.path.isfile(os.path.join(source, item)):
                self.put(os.path.join(source, item), '%s/%s' % (target, item))
            else:
                self.mkdir('%s/%s' % (target, item), ignore_existing=True)
                self.put_dir(os.path.join(source, item), '%s/%s' % (target, item))
    def mkdir(self, path, mode=751, ignore_existing=False):
        ''' Augments mkdir by adding an option to not fail if the folder exists  '''
        try:
            super(MySFTPClient, self).mkdir(path, mode)
        except IOError:
            if ignore_existing:
                pass
            else:
                raise

pipelines = ["ATAC-seq","ChIP-seq","Cut_Tag","RNA-seq","WGBS"]

clusters = {
    "hpc2.ust.hk" : "/home/share/dcyleung/snakemake/", 
    "biocrfhpc1.ust.hk" : "/home4/share/dcyleung/snakemake/", 
    "biocrfhpc2.ust.hk" : "/data1/share/dcyleung/Pipeline/snakemake/"
    }

parser = argparse.ArgumentParser(
    description="Updates specified pipelines in all cluster",
)

parser.add_argument(
    '--pipes', '-p',
    type=str,
    nargs='+',
    required=True,
    help = "Pipelines to update separated by whitespaces. Has to be exact names"
    )

args = parser.parse_args()

if "all" in args.pipes:
    args.pipes = pipelines

for pipe in args.pipes:
    if pipe not in pipelines:
        sys.exit(f"The pipe {pipe} is not a possible choice to update the pipelines")

for pipe in args.pipes:
    for cluster in clusters.keys():
        transport = paramiko.Transport((cluster))
        transport.connect(username=os.environ['SSH_HKUST_USER'],password=os.environ['SSH_HKUST_PASS'])
        sftp = MySFTPClient.from_transport(transport)
        source_path = pipe + "/workflow/"
        target_path = clusters[cluster] + source_path
        #For the folder you refuse to change the name
        if source_path == "ChIP-seq/workflow/":
            source_path = "Chipseq/workflow/"
        sftp.put_dir(source_path, target_path)
        sftp.close()
        print(f"Updated pipeline {pipe} in cluster {cluster}")