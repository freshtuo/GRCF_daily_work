#!/usr/bin/env python
# send_notification_email
# send a notification email to users for their data
# Version 1.0
# Author: Tuo Zhang
# Date: 06/23/2021
# NEW: 
# 

import sys
import smtplib
import mimetypes
import yaml

from os.path import exists
from os.path import basename
from time import time
from time import localtime
from time import asctime
from email.message import EmailMessage
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter

# class
class MyEmail:
    """Class for sending notification email"""

    def __init__(self, args):
        """class constructor"""
        # dictionary saving email settings
        self.setdic = {}
        # email main text
        self.main_text = ''
        # complete email message
        self.msg = EmailMessage()
        # commandline arguments
        self.args = args

    def get_time(self):
        """get current time"""
        return '['+asctime(localtime(time()))+']'

    def initialize_settings(self):
        """load settings from a yaml file"""
        setting_file = self.args.setting
        # read from setting yaml if exists
        if exists(setting_file):
            with open(setting_file, 'r') as fset:
                self.setdic = yaml.safe_load(fset)
        else:
            # initialize settings
            # 'run' section: sequencing run info
            self.setdic['run'] = {'iLab':-1, 'libType':'', 'nSamples':0, 'instrument':'', 'seqType':'', 'readLen':[], 'date':''}
            # 'user' section: user contact info
            self.setdic['user'] = {'name':'', 'email':'', 'piEmail':'', 'dataPath':''}
            # 'core' section: core contact info
            self.setdic['core'] = {'name':'Adrian', 'fromEmail':'yit2001@med.cornell.edu', 'ccEmails':['yit2001@med.cornell.edu', 'taz2008@med.cornell.edu']}
            # 'attachments' section: files to attach
            self.setdic['attachments'] = {'demuxSum':[], 'other':[]}

    def prepare_main_text(self):
        """prepare email main text"""
        # collect info
        instrument = self.setdic['run']['instrument']
        ilab = self.setdic['run']['iLab']
        nsamples = self.setdic['run']['nSamples']
        dataloc = self.setdic['user']['dataPath']
        username = self.setdic['user']['name']
        sender = self.setdic['core']['name']
        # update main text
        if self.args.location == 'sftp':
            self.main_text = """\
<html>
  <head></head>
  <body>
    <p>Hi {},</p>
    <p>Your {} sequencing run ({} samples, iLab request: {}) is complete. 
       The data have been uploaded to your lab folder on our internal data transferring server, 
       which enables both ssh and sftp access with your Cornell CWID.</p>
    <p>Here are the account information:</p>
    <p>ssh(sftp) CWID@157.139.241.40<br>
       Account: CWID<br>
       Password: password associated with CWID</p>
    <p>The data can be downloaded from the following directories:</p>
    <p style="color:blue;font-size:18px;"><b>{}</b></p>
    <p>You can certainly use terminal to access it but a ftp client like Filezilla (https://filezilla-project.org/) might make your life easier.</p>
    <p><b>The data will be available for 4 weeks, after that they will be removed from the server. 
          Please download your data on time.</b></p>
    <p>The run summary was attached. Please let us know if there are any problems.</p>
    <p>Best,</p>
    <p>{}</p>
  </body>
</html>
""".format(username,instrument,nsamples,ilab,dataloc,sender)
        elif self.args.location == 'aws':
            self.main_text = """\
<html>
  <head></head>
  <body>
    <p>Hi {},</p>
    <p>Your {} sequencing run ({} samples, iLab request: {}) is complete.
       The data can be downloaded via the following AWS s3 URL:</p>
    <p style="color:blue;font-size:18px;"><b>{}</b></p>
    <p>You can either download the data with a regular web browser by clicking the above link 
       or through a Linux terminal with a command like this one: </p>
    <p>wget -O YOU-PREFERED-FILE-NAME.tar  “THE-ABOVE-URL"</p>
    <p><b>The data will be available for 7 days only, after that they will be removed. Please download your data on time.</b></p>
    <p>The run summary was attached. Please let us know if there are any problems.</p>
    <p>Best,</p>
    <p>{}</p>
  </body>
</html>
""".format(username,instrument,nsamples,ilab,dataloc,sender)

    def prepare_email(self):
        # collect info
        instrument = self.setdic['run']['instrument']
        date = self.setdic['run']['date']
        libtype = self.setdic['run']['libType']
        seqtype = self.setdic['run']['seqType']
        readlen = '+'.join(['{}'.format(x) for x in self.setdic['run']['readLen']])
        e_tos = [self.setdic['user']['email'], self.setdic['user']['piEmail']]
        e_ccs = self.setdic['core']['ccEmails']
        e_from = self.setdic['core']['fromEmail']
        # demultiplex summary files
        demuxsum = []
        if 'demuxsum' in self.setdic['attachments']:
            for sum_file in self.setdic['attachments']['demuxSum']:
             # file exists?
             if exists(sum_file):
                 demuxsum.append(sum_file)
        # other files
        other = []
        if 'other' in self.setdic['attachments']:
            for other_file in self.setdic['attachments']['other']:
             # file exists?
             if exists(other_file):
                 other.append(other_file)
        # create email message
        self.msg['Subject'] = '{} {}{} {} {} sequencing data'.format(instrument,seqtype,readlen,date,libtype)
        self.msg['To'] = ', '.join(e_tos)
        self.msg['Cc'] = ', '.join(e_ccs)
        self.msg['From'] = e_from
        # email content
        self.msg.set_content(self.main_text, subtype='html')
        # add attachments (demux summary & other files)
        for dfile in demuxsum+other:
            ctype, encoding = mimetypes.guess_type(dfile)
            if ctype is None or encoding is not None:
                # No guess could be made, or the file is encoded (compressed),
                # so use a generic bag-of-bits type.
                ctype = 'application/octet-stream'
            maintype, subtype = ctype.split('/', 1)
            with open(dfile, 'rb') as fp:
                self.msg.add_attachment(fp.read(), maintype=maintype, subtype=subtype, filename=basename(dfile))

    def send_email():
        """"send the message via local SMTP server."""
        with smtplib.SMTP('localhost') as s:
            s.send_message(self.msg)

# functions
def get_arguments():
    """fetch commandline arguments."""
    parser = ArgumentParser(description="""Send a notification email to user""", formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-s', '--setting', required=True, help="""a yaml file specifying detailed information.""")
    parser.add_argument('-l', '--location', choices=['sftp','aws'], default='sftp', help="""data location either sftp or aws.""")
    parser.add_argument('-f', '--fastq', help="""path to the fastq files""")
    return parser.parse_args()

def main():
    """call me to get started!"""
    args = get_arguments()
    m = MyEmail(args)
    m.initialize_settings()
    m.prepare_main_text()
    m.prepare_email()
    print(m.msg)
    #m.send_email()

# main
if __name__ == '__main__':
    main()
    #with open('/tmp/1.yaml', 'w') as fout:
    #    yaml.dump(setdic, fout, default_flow_style=False)

