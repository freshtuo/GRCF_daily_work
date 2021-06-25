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
# get current time
def get_time():
	return '['+asctime(localtime(time()))+']'

# extract information out of the setting yaml file
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

# prepare email main text based on data location
def prepare_content(setdic, location):
	content = ''
	if location == 'sftp':
		# collect info
		instrument = setdic['run']['instrument']
		ilab = setdic['run']['iLab']
		nsamples = setdic['run']['nSamples']
		dataloc = setdic['user']['dataPath']
		username = setdic['user']['name']
		sender = setdic['core']['name']
		# update text
		content = """\
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
	elif location == 'aws':
		# collect info
		instrument = setdic['run']['instrument']
		ilab = setdic['run']['iLab']
		nsamples = setdic['run']['nSamples']
		dataloc = setdic['user']['dataPath']
		username = setdic['user']['name']
		sender = setdic['core']['name']
		# update text
		content = """\
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
	return content

def prepare_email(setting_file, location):
	# parse settings
	setdic = parse_settings(setting_file)
	# collect info
	instrument = setdic['run']['instrument']
	date = setdic['run']['date']
	libtype = setdic['run']['libType']
	seqtype = setdic['run']['seqType']
	readlen = '+'.join(['{}'.format(x) for x in setdic['run']['readLen']])
	demuxsum = []
	if setdic['attachments'] is not None:
		demuxsum = setdic['attachments']['demuxSum']
	e_tos = [setdic['user']['email'], setdic['user']['piEmail']]
	e_ccs = setdic['core']['ccEmails']
	e_from = setdic['core']['fromEmail']
	# create email message
	msg = EmailMessage()
	msg['Subject'] = '{} {}{} {} {} sequencing data'.format(instrument,seqtype,readlen,date,libtype)
	msg['To'] = ', '.join(e_tos)
	msg['Cc'] = ', '.join(e_ccs)
	msg['From'] = e_from
	# email content
	content = prepare_content(setdic, location)
	msg.set_content(content, subtype='html')
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

def send_email(setting_file, location):
	msg = prepare_email(setting_file, location)
	# Send the message via local SMTP server.
	with smtplib.SMTP('localhost') as s:
		s.send_message(msg)

def get_arguments():
	parser = ArgumentParser(description="""Send a notification email to user""")
	parser.add_argument('-s', '--setting', required=True, help="""a yaml file specifying detailed information.""")
	parser.add_argument('-l', '--location', choices=['sftp','aws'], default='sftp', help="""data location either sftp or aws.""")
	return parser.parse_args()

if __name__ == '__main__':
	args = get_arguments()
	send_email(args.setting, args.location)

