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

import xml.etree.ElementTree as ET

from os import listdir
from os.path import exists
from os.path import basename
from os.path import isdir
from time import time
from time import localtime
from time import asctime
from email.message import EmailMessage
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter
from re import search

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
        # available sequencers
        self.instruments = ['NovaSeq6000','HiSeq4000','NextSeq500','NextSeq2000','MiSeq']
        # readlength tags used in xml file
        self.tagdic = {'NovaSeq6000':('Read1NumberOfCycles','Read2NumberOfCycles'),
            'HiSeq4000':('Read1','Read2'),
            'NextSeq500':('Read1','Read2'),
            'NextSeq2000':('Read1','Read2'),
            'MiSeq':('RunInfoRead','NumCycles')}
        # fastq path (in order to remove the last '/' character if any)
        self.fastq_path = ''

    def get_time(self):
        """get current time"""
        return '['+asctime(localtime(time()))+']'

    def initialize_settings(self):
        """load settings from a yaml file"""
        setting_file = self.args.setting
        # read from setting yaml if exists
        if exists(setting_file):
            print('load settings from file {}\n'.format(setting_file))
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
        # user assigned library type?
        if self.args.library is not None:
            self.setdic['run']['libType'] = self.args.library
        # user assigned fastq path?
        if self.args.fastq is not None:
            self.fastq_path = self.args.fastq
            # remove the last '/' if it appears at the end of the fastq path
            if self.fastq_path[-1] == '/':
                self.fastq_path = self.fastq_path[:-1]
            

    def infer_instrument(self):
        """guess which sequencer was used"""
        for ins in self.instruments:
            if ins in self.fastq_path:
                self.setdic['run']['instrument'] = ins
                print('Infer instrument: {}'.format(ins))
                break

    def infer_date(self):
        """guess sequencing date"""
        # e.g. '210615_A00814_0436_BH7VLVDSX2'
        tpat = search('/(\d+)_[a-zA-Z0-9\-]+_\d+', self.fastq_path)
        if tpat:
            date = tpat.groups()[0]
            self.setdic['run']['date'] = '20{}_{}_{}'.format(date[:2], date[2:4], date[4:])
            print('Infer date: {}'.format(self.setdic['run']['date']))

    def infer_ilab(self):
        """guess iLab id"""
        # Loda-MJ-10557_2021_06_15
        items = self.fastq_path.split('/')[-1].split('-')
        if len(items) > 2:
            tpat = search('(\d+)',items[2])
            if tpat:
                ilab = int(tpat.groups()[0])
                self.setdic['run']['iLab'] = ilab
                print('Infer iLab id: {}'.format(ilab))

    def infer_run_folder(self):
        """guess sequencing run folder"""
        # locate 'Unaligned' in fastq path
        items = self.fastq_path.split('/')
        pos = -1
        for k,x in enumerate(items):
            if 'Unaligned' in x:
                pos = k
                break
        # found!
        if pos != -1:
            return '/'.join(items[:pos])
        else:
            return None

    def infer_run_para_file(self):
        """guess RunParameters.xml file location"""
        run_folder = self.infer_run_folder()
        if run_folder is not None:
            for tfile in listdir(run_folder):
                if tfile.upper() == 'RunParameters.xml'.upper():
                    return '{}/{}'.format(run_folder,tfile)
        return None

    def infer_seqinfo(self):
        """guess sequencing type and read length"""
        # tag available for the given sequencer?
        instrument = self.setdic['run']['instrument']
        if instrument not in self.tagdic:
            return None
        # fetch tag name for read 1 & 2
        tag1,tag2 = self.tagdic[instrument]
        # try to find RunParameters.xml file
        run_para_file = self.infer_run_para_file()
        ##print(run_para_file)
        if run_para_file is not None:
            # parse xml file
            tree = ET.parse(run_para_file)
            # get root
            root = tree.getroot()
            # fetch read length
            read_length = []
            if instrument == 'MiSeq':
                # e.g.
                # <RunInfoRead Number="1" NumCycles="151" IsIndexedRead="N" />
                # <RunInfoRead Number="2" NumCycles="8" IsIndexedRead="Y" />
                for x in root.iter('RunInfoRead'):
                    if x.attrib['IsIndexedRead'] == 'N':# sample read, not index!
                        if x.attrib['NumCycles'] != '0':
                            read_length.append(int(x.attrib['NumCycles']))
            else:# not MiSeq
                read_length = []
                # get read 1 length
                for x in root.iter(tag1):
                    if x.text != '0':
                        read_length.append(int(x.text))
                    break
                # get read 2 length
                for x in root.iter(tag2):
                    if x.text != '0':
                        read_length.append(int(x.text))
                    break
            # update info in settings
            self.setdic['run']['readLen'] = read_length
            if len(read_length) > 1:
                self.setdic['run']['seqType'] = 'PE'
            else:
                self.setdic['run']['seqType'] = 'SR'
            print('infer seqType: {}'.format(self.setdic['run']['seqType']))
            print('infer readLen: {}'.format(read_length))
        return None

    def infer_nsamples(self):
        """guess number of sequenced samples"""
        # search for fastq files
        # e.g. 10_S10_L001_R1_001.fastq.gz
        sids = list(set([search('.*_S(\d+)_L00\d+_R1_001',tfile).groups()[0] for tfile in listdir(self.fastq_path) if search('R1_001\.fastq\.gz',tfile)]))
        # search for folders other than 'Summary' if failing to find *.fastq.gz
        if not sids:
            sids = [tfile for tfile in listdir(self.fastq_path) if isdir('{}/{}'.format(self.fastq_path,tfile)) and tfile != 'Summary']
        # update info in settings
        if sids:
            self.setdic['run']['nSamples'] = len(sids)
            print('infer nSamples: {}'.format(len(sids)))

    def infer_demuxsum(self):
        """guess demux summary html files"""
        # get flowcell id
        run_folder = self.infer_run_folder()
        fcid = run_folder.split('/')[-1].split('_')[-1]
        if self.setdic['run']['instrument'] not in ['NextSeq2000','MiSeq']:
            fcid = fcid[1:]# the first letter refers to flowcell A or B
        # get user folder
        user_folder = self.fastq_path.split('/')[-1]
        # locate report folder
        demuxs = []
        for x in listdir('{}/Reports/html/{}'.format('/'.join(self.fastq_path.split('/')[:-1]), fcid)):
            if x in user_folder:# find it!
                demuxs.append('{}/Reports/html/{}/{}/all/all/lane.html'.format('/'.join(self.fastq_path.split('/')[:-1]), fcid, x))
                demuxs.append('{}/Reports/html/{}/{}/all/all/laneBarcode.html'.format('/'.join(self.fastq_path.split('/')[:-1]), fcid, x))
        # update info in settings
        if demuxs:
            self.setdic['attachments']['demuxSum'] = demuxs
            print('infer demux summary: {}'.format(demuxs))

    def infer_settings(self):
        """infer settings based on fastq path, and overwrite the current one"""
        # fastq location assigned by user?
        if self.args.fastq is None:
            return None
        # fastq location exists?
        if not exists(self.fastq_path):
            return None
        print('infer settings based on fastq path {}'.format(self.args.fastq))
        # instrument
        self.infer_instrument()
        # date
        self.infer_date()
        # iLab
        self.infer_ilab()
        # seqType & readLen
        self.infer_seqinfo()
        # nSamples
        self.infer_nsamples()
        # demuxSum
        self.infer_demuxsum()
        print('')

    def write_settings(self):
        """write settings to file"""
        setting_file = self.args.setting
        # write only if:
        # 1) the setting file does not exist
        # 2) or --overwrite option is on
        if exists(setting_file) and (not self.args.overwrite):
            print('setting file already exists. use --overwrite to force writing to it!\n')
            return None
        else:
            if not exists(setting_file):
                print('write settings to file {}\n'.format(setting_file))
            else:
                print('overwrite setting file {}\n'.format(setting_file))
            with open(setting_file, 'w') as fout:
                yaml.dump(self.setdic, fout, default_flow_style=False)

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
        e_tos = list(set([self.setdic['user']['email'], self.setdic['user']['piEmail']]))
        e_ccs = self.setdic['core']['ccEmails']
        e_from = self.setdic['core']['fromEmail']
        # demultiplex summary files
        demuxsum = []
        if 'demuxSum' in self.setdic['attachments']:
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
        if self.args.note is None:
            self.msg['Subject'] = '{} {}{} {} {} sequencing data'.format(instrument,seqtype,readlen,date,libtype)
        else:
            self.msg['Subject'] = '{} {}{} {} {} sequencing data {}'.format(instrument,seqtype,readlen,date,libtype,self.args.note)
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

    def check_settings(self):
        """check settings before sending out email."""
        myPass = True
        # run section
        if self.setdic['run']['iLab'] == -1:
            print('warning: iLab id is missing.')
        if self.setdic['run']['libType'] == '':
            print('warning: library type is missing.')
        if self.setdic['run']['nSamples'] == 0:
            print('warning: number of samples is missing.')
        if self.setdic['run']['instrument'] == '':
            print('warning: sequencer info is missing.')
        if not self.setdic['run']['readLen']:
            print('warning: read length is missing.')
        if self.setdic['run']['date'] == '':
            print('warning: sequencing date is missing.')
        # user section
        for item in ['name','email','piEmail','dataPath']:
            if not self.setdic['user'][item]:
                print('error: user {} is missing.'.format(item))
                myPass = False
        # core section (skipped)
        # attachments section
        if not self.setdic['attachments']['demuxSum']:
            print('warning: demux summary files are missing.')
        print('')
        return myPass

    def print_main_text(self):
        """show email main text."""
        if self.args.text:
            print('##################################### main text ######################################')
            print(self.main_text)
            print('######################################################################################\n')

    def print_msg(self):
        """print full email message."""
        if self.args.msg:
            print('################################### full email msg ###################################')
            print(self.msg)
            print('######################################################################################\n')

    def __repr__(self):
        """print settings."""
        print('###################################### settings ######################################')
        for section in ['run','user','core','attachments']:
            print('{}:'.format(section))
            for item in self.setdic[section]:
                print('\t{}: {}'.format(item, self.setdic[section][item]))
        print('######################################################################################\n')

    def print_settings(self):
        """print settings."""
        if self.args.print:
            self.__repr__()

    def send_email(self):
        """send the message via local SMTP server."""
        # send out email only if passing settings check and '--email' option is turned on
        if self.args.email and self.check_settings():
            with smtplib.SMTP('localhost') as s:
                s.send_message(self.msg)
        else:
            print('once ready, please use --email option to send this email out!\n')

# functions
def get_arguments():
    """fetch commandline arguments."""
    parser = ArgumentParser(description="""Send a notification email to user""", formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-s', '--setting', required=True, help="""a yaml file specifying detailed information.""")
    parser.add_argument('-l', '--location', choices=['sftp','aws'], default='sftp', help="""data location either sftp or aws.""")
    parser.add_argument('-b', '--library', default='RNAseq', help="""library type.""")
    parser.add_argument('-f', '--fastq', help="""absolute path to the fastq files, do not put '/' at the end!!!""")
    parser.add_argument('-w', '--overwrite', action='store_true', help="""allow overwriting setting file.""")
    parser.add_argument('-e', '--email', action='store_true', help="""allow sending out email.""")
    parser.add_argument('-t', '--text', action='store_true', help="""print email main text.""")
    parser.add_argument('-m', '--msg', action='store_true', help="""print full email msg.""")
    parser.add_argument('-p', '--print', action='store_true', help="""print settings.""")
    parser.add_argument('-n', '--note', help="""include a note at the end of email subject.""")
    return parser.parse_args()


def main():
    """call me to get started!"""
    args = get_arguments()
    m = MyEmail(args)
    m.initialize_settings()
    m.infer_settings()
    m.write_settings()
    m.prepare_main_text()
    m.prepare_email()
    ##print(m.check_settings())
    m.print_main_text()
    m.print_msg()
    m.print_settings()
    m.send_email()

# main
if __name__ == '__main__':
    main()

