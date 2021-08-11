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

