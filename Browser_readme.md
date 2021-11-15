# Index
- [Index](#index)
- [Introduction](#introduction)
- [Transferring the files](#transferring-the-files)
  - [sftp](#sftp)
  - [scp](#scp)
- [Loading the bw files in UCSC genome browser](#loading-the-bw-files-in-ucsc-genome-browser)

# Introduction

This readme provides instructions on how to use the `.bw` files obtained as a result of a pipeline to visualize your results in the genome browser.

To use it, you first need an account in the lab synology server, and know where your `.bw` files are located in the cluster.

# Transferring the files

First login in the cluster and go to the folder where the `.bw` files are located, then you have 2 options to transfer the files, `sftp` or `scp`.

## sftp

sftp allows you to connect remotely and open a session in the synology server, which allows you to make folders, move files, etc. like you would in the cluster but in synology. To connect to synology you can use:

    sftp <username>@lfz097.ust.hk

**REPLACING \<username> BY YOUR OWN USERNAME** in synology and using your synology password. Once you are connected to synology, you can create or modify the directory where you want to store your `.bw` files. Once you are happy with your location you can use

    put <your bw file>

to transfer a bw file to the synology from the cluster, replacing the text within \<> with your file name. Remember that you can also use wildcards (*), so if you want to transfer all `.bw` files in the cluster directory to the synology you can use:

    put *.bw

Once you are done transferring or modifying files on the synology, you can use

    exit

to leave the synolgy session.

## scp

Unlike sftp, scp only transfers files between the server and synology, so it does not allow you to create or modify files and folders in synology. To transfer files using scp you can do:

    scp /path/to/files/in/cluster <username>@lfz097.ust.hk:~/your/synology/folder/

As before, you can use wildcards, so to transfer all `.bw` files in a directory you can use for example:

    scp /path/to/files/in/cluster/*.bw <username>@lfz097.ust.hk:~/your/synology/folder/

The files should then be transferred to the synology.

# Loading the bw files in UCSC genome browser

Once your files are located in the synology server, they can be accessed online, so we can use the UCSC genome browser to website to upload your data into the genome browser. First go to the UCSC website, and then to the custom tracks section or just this url:

> https://genome.ucsc.edu/cgi-bin/hgCustom

First, select the correct genome for UCSC genome browser at the top of the page, and then you can paste the `bw` file url in the "Paste URL or data" dialogue box. You can use the following template to upload your files:

    track type=bigWig name="Your_sample_name" description="Your_sample_description" color=255,211,0 visibility=2 windowingFunction=maximum autoScale=on maxHeightPixels=30 bigDataUrl=https://lfz097.ust.hk/~<username>/your/bw/folder/sample.bw

and click on "submit" to upload your data. Be mindful of new lines or misspelling mistakes in this line, as they will cause the upload to fail and UCSC will complain.