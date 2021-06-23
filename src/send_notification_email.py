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

# functions
def get_time():
	return '['+asctime(localtime(time()))+']'

def parse_settings(setting_file):
	# read in setting info
	with open(setting_file, 'r') as fset:
		setdic = yaml.safe_load(fset)
	# check existence of attachments
	if setdic['attachments'] == None:
		print('No attachments detected')
	else:
		for sum_file in setdic['attachments']['demuxSum']:
			if not exists(sum_file):
				print('Error: Cannot locate demux summary file {}!'.format(sum_file))
				sys.exit(1)
	return(setdic) 

def prepare_email(setting_file):
	# parse settings
	setdic = parse_settings(setting_file)
	# collect info
	instrument = setdic['run']['instrument']
	date = setdic['run']['date']
	libtype = setdic['run']['libType']
	ilab = setdic['run']['iLab']
	nsamples = setdic['run']['nSamples']
	seqtype = setdic['run']['seqType']
	readlen = '+'.join(['{}'.format(x) for x in setdic['run']['readLen']])
	demuxsum = []
	if setdic['attachments'] is not None:
		demuxsum = setdic['attachments']['demuxSum']
	dataloc = setdic['user']['dataPath']
	username = setdic['user']['name']
	e_tos = [setdic['user']['email'], setdic['user']['piEmail']]
	e_ccs = setdic['core']['ccEmails']
	e_from = setdic['core']['fromEmail']
	sender = setdic['core']['name']
	# create email message
	msg = EmailMessage()
	msg['Subject'] = '{} {}{} {} {} sequencing data'.format(instrument,seqtype,readlen,date,libtype)
	msg['To'] = ', '.join(e_tos)
	msg['Cc'] = ', '.join(e_ccs)
	msg['From'] = e_from
	# email content
	msg.set_content("""\
Hi {},

Your {} sequencing run ({} samples, iLab request: {}) is complete. \
The data have been uploaded to your lab folder on our internal data transferring server, \
which enables both ssh and sftp access with your Cornell CWID.

Here are the account information: 
ssh(sftp) CWID@157.139.241.40 
Account: CWID 
Password: password associated with CWID 
The data can be downloaded from the following directories: 

{}

You can certainly use terminal to access it but a ftp client like Filezilla (https://filezilla-project.org/) might make your life easier. 

The data will be kept there for about one month. 

The run summary was attached. Please let us know if there are any problems. 
	 
Best, 
{}
""".format(username,instrument,nsamples,ilab,dataloc,sender))
	# add demux summary if exists
	for dfile in demuxsum:
		ctype, encoding = mimetypes.guess_type(dfile)
		if ctype is None or encoding is not None:
			# No guess could be made, or the file is encoded (compressed), so
			# use a generic bag-of-bits type.
			ctype = 'application/octet-stream'
		maintype, subtype = ctype.split('/', 1)
		with open(dfile, 'rb') as fp:
			msg.add_attachment(fp.read(), maintype=maintype, subtype=subtype, filename=basename(dfile))
	return msg

def send_email(setting_file):
	msg = prepare_email(setting_file)
	# Send the message via local SMTP server.
	with smtplib.SMTP('localhost') as s:
		s.send_message(msg)

def get_arguments():
	parser = ArgumentParser(description="""Send a notification email to user""")
	parser.add_argument('-s', '--setting', required=True, help="""a yaml file specifying detailed information.""")
	return parser.parse_args()

if __name__ == '__main__':
	args = get_arguments()
	send_email(args.setting)

