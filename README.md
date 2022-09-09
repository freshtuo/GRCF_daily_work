# GRCF_daily_work

scripts in hope of simplifying daily work at GRCF

1. script for preparing and sending out notification emails to users

<pre>
usage: send_notification_email.py [-h] -s SETTING [-l {sftp,aws}] [-b LIBRARY]
                                  [-f FASTQ] [-w] [-e] [-t] [-m] [-p]
                                  [-n NOTE]

Send a notification email to user

optional arguments:
  -h, --help            show this help message and exit
  -s SETTING, --setting SETTING
                        a yaml file specifying detailed information. (default:
                        None)
  -l {sftp,aws}, --location {sftp,aws}
                        data location either sftp or aws. (default: sftp)
  -b LIBRARY, --library LIBRARY
                        library type. (default: None)
  -f FASTQ, --fastq FASTQ
                        absolute path to the fastq files, do not put '/' at
                        the end!!! (default: None)
  -w, --overwrite       allow overwriting setting file. (default: False)
  -e, --email           allow sending out email. (default: False)
  -t, --text            print email main text. (default: False)
  -m, --msg             print full email msg. (default: False)
  -p, --print           print settings. (default: False)
  -x SUFFIX, --suffix SUFFIX
                        add a suffix to the end of email subject (default:
                        None)
  -n NOTE, --note NOTE  include a note to the email main text. (default: )
</pre>

2. script for examining/spliting samplesheet and generating shell scripts for
   demultiplexing with bcl2fastq2

<pre>
usage: auto_bcl2fq_script.py [-h] [-v] [-r [run_folder]] -i [samplesheet_file]
                             [-s [lanes_to_skip [lanes_to_skip ...]]]
                             [-b [bcl2fastq]] [-o [script_folder]]
                             [-e [samplesheet_folder]] [-f] [-n]
                             [-p [{hiseq,nextseq500,nextseq2000,novaseq}]]

Given a samplesheet, check its format and automatically generate shell script
for demux.

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -r [run_folder], --runfolder [run_folder]
                        sequencing run folder
  -i [samplesheet_file], --samplesheet [samplesheet_file]
                        samplesheet file
  -s [lanes_to_skip [lanes_to_skip ...]], --skiplanes [lanes_to_skip [lanes_to_skip ...]]
                        skip lanes (e.g. -s 2 3 4)
  -b [bcl2fastq], --bcl2fastq [bcl2fastq]
                        location of bcl2fastq program
  -o [script_folder], --scriptdir [script_folder]
                        output folder to write shell script(s) (default: run_folder)
  -e [samplesheet_folder], --samplesheetdir [samplesheet_folder]
                        output folder to write samplesheet(s) (default: basecall_folder)
  -f, --force           whether or not to overwrite output script if existing
  -n, --no-lane-splitting
                        whether or not to add --no-lane-splitting option to shell script
  -p [{hiseq,nextseq500,nextseq2000,novaseq}], --platform [{hiseq,nextseq500,nextseq2000,novaseq}]
                        sequencing platform
</pre>

3. script to examining conflicting index sequences in the samplesheet file
   before mixing samples and loading them to the instrument

<pre>
usage: review_samplesheet_index.py [-h] [-v] -i [samplesheet_file] [-l
[log_file]] [-d] [-o [updated_samplesheet_file]] [-t [tenx_dir]] [-r
[tenx_idx_file]]

Given a samplesheet, check index sequence conflicts.

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -i [samplesheet_file], --samplesheet [samplesheet_file]
                        samplesheet file
  -l [log_file], --log [log_file]
                        log file
  -d, --detail          output details
  -o [updated_samplesheet_file], --updatedsheet [updated_samplesheet_file]
                        updated samplesheet file for manual check
  -t [tenx_dir], --tenxdir [tenx_dir]
                        load 10X index sequences in the given folder
  -r [tenx_idx_file], --tenxidx [tenx_idx_file]
                        load 10X index sequences from the given file
</pre>

