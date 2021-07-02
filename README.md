# GRCF_daily_work

scripts in hope of simplifying daily work at GRCF

1. script for preparing and sending out notification emails to users

usage: send_notification_email.py [-h] -s SETTING [-l {sftp,aws}] [-b LIBRARY]\
&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; [-f FASTQ] [-w] [-e] [-t] [-m] [-p]\
&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; [-n NOTE]

Send a notification email to user

optional arguments:\
&nbsp; -h, --help&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; show this help message and exit\
&nbsp; -s SETTING, --setting SETTING\
&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; a yaml file specifying detailed information. (default:
&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; None)
&nbsp; -l {sftp,aws}, --location {sftp,aws}\
&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; data location either sftp or aws. (default: sftp)\
&nbsp; -b LIBRARY, --library LIBRARY\
&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; library type. (default: None)\
&nbsp; -f FASTQ, --fastq FASTQ\
&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; absolute path to the fastq files, do not put '/' at
&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; the end!!! (default: None)
&nbsp; -w, --overwrite&nbsp; &nbsp; &nbsp;  allow overwriting setting file. (default: False)\
&nbsp; -e, --email&nbsp; &nbsp; &nbsp; &nbsp; &nbsp;  allow sending out email. (default: False)\
&nbsp; -t, --text&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; print email main text. (default: False)\
&nbsp; -m, --msg&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;  print full email msg. (default: False)\
&nbsp; -p, --print&nbsp; &nbsp; &nbsp; &nbsp; &nbsp;  print settings. (default: False)\
&nbsp; -n NOTE, --note NOTE&nbsp; include a note at the end of email subject. (default:
&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; None)

