#!/usr/bin/env python
# arrange_fastq.py
# arrange demux fastq files for distribution
# 

import sys
import re
import pandas as pd
import os
import logging
import subprocess

from collections import Counter
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter

# class
class MyBCLConvert:
	"""Class for handling demux data processed by Illumina BCL Convert"""
	def __init__(self, args):
		"""class constructor"""
		# samplesheet filename
		self.sample_sheet_file = args.samplesheetfile
		# info sections in samplesheet
		self.sections = {}
		# number of rows in samplesheet
		self.nrow = None
		# table that stores sample by lane (DataFrame)
		self.lane_info = None
		# table that stores sample by project (DataFrame)
		self.project_info = None
		# samplesheet in v1 format (DataFrame)
		self.v1_format = None
		# v1 format samplesheet file
		self.v1_format_file = args.v1formatfile
		# demux folder
		self.demux_folder = args.demuxfolder
		# demux stats file
		self.demux_stats_file = None
		# demux BCLConvert folder
		self.demux_bclconvert_folder = None
		# demux fastq folder
		self.demux_fastq_folder = None
		# demux stats (DataFrame)
		self.demux_stats = None
		# demux fastq files (DataFrame)
		self.demux_fastq_files = None
		# demux fastqc report files (DataFrame)
		self.demux_fastqc_report_files = None
		# output folder
		self.output_folder = args.outfolder
		# demux stats output columns
		self.demux_cols = ['Lane','SampleID','Index','# Reads','% Perfect Index Reads','% One Mismatch Index Reads','ProjectName']
		# v1 format samplesheet output columns
		self.v1_format_cols = ['Lane','Sample_ID','Sample_Name','Sample_Plate','Sample_Well','Index_Plate_Well',\
		'I7_Index_ID','index','I5_Index_ID','index2','Sample_Project','Description']
		# shell script for arranging fastq files
		self.arrange_fastq_script = None
		# shell script for arranging fastqc report files
		self.arrange_fastqc_report_script = None
		# shell script for running multiqc to merge fastqc results from multiple samples
		self.run_multiqc_script = None
		self.run_multiqc_script_log = None
		# make symbolic links instead of hard links
		self.symbolic = args.symbolic
		# generate link-file shell scripts without running them
		self.norun = args.norun
		# overall demux summary (DataFrame)
		self.demux_summary = None
		# overall demux summary output file
		self.demux_summary_file = args.sumfile
		# multiqc executable
		self.multiqc = args.multiqc
		# run read cycle info (Dict)
		self.cycles = Counter()
		# run read cycle labels
		self.read_1_label = 'Read1Cycles'
		self.read_2_label = 'Read2Cycles'
		self.index_1_label = 'Index1Cycles'
		self.index_2_label = 'Index2Cycles'
		# run override cycles (in case not provided in sample sheet)
		self.override_cycles = None
		# skip projects for multiqc
		self.skip_multqc = args.skipmultqc

	def split_sections(self):
		"""preprocess sample sheet: split info by sections"""
		last_pos = -1
		section_name = ''
		with open(self.sample_sheet_file, 'r') as fin:
			data = fin.readlines()
			for k,line in enumerate(data):
				tpat = re.search("^\[(.*?)\]", line)
				if tpat:
					if last_pos != -1:
						self.sections[section_name] = (last_pos, k)
						last_pos = k
					else:
						last_pos = k
					section_name = tpat.groups()[0]
			# add the last section
			self.sections[section_name] = (last_pos, len(data))
			# record the total number of rows in sample sheet
			self.nrow = len(data)
			# logging
			logging.info('{} sections detected.'.format(len(self.sections)))

	def infer_read_cycles(self):
		"""infer read cycles from [Reads] section"""
		if 'Reads' not in self.sections:
			logging.error('Section Reads is unavailable.')
		start, end = self.sections['Reads']
		with open(self.sample_sheet_file, 'r') as fin:
			for entry in fin.readlines()[start+1:end-1]:
				attr = entry.split(',')[0].strip()
				value = int(entry.split(',')[1].strip())
				self.cycles[attr] = value
		# self.cycles Counter(), thus no need to assign zero for missing read/index labels
		self.override_cycles = 'Y{};I{};I{};Y{}'.format(self.cycles[self.read_1_label], self.cycles[self.index_1_label], \
			self.cycles[self.index_2_label], self.cycles[self.read_2_label])
		logging.info('Infer run read cycles: {}.'.format(self.override_cycles))

	def load_section_table(self, section_name):
		"""load data from a given section into a DataFrame"""
		if section_name not in self.sections:
			logging.error('Section {} is unavailable.'.format(section_name))
			sys.exit(2)
		# read in table
		start, end = self.sections[section_name]
		df = pd.read_table(self.sample_sheet_file, sep=',', header=0, skiprows=start+1, skipfooter=self.nrow-end, engine='python')
		# remove empty rows
		df.dropna(axis=0, how='all', inplace=True)
		# logging
		logging.info('Load section {}: {}.'.format(section_name, df.shape))
		return df

	def get_sample_info(self):
		"""extract project and sample information"""
		logging.info('Extracting project and sample information.')
		# file exists?
		if not os.path.exists(self.sample_sheet_file):
			logging.error('failed.\nUnable to locate the samplesheet file: {}.'.format(self.samplesheetfile))
			sys.exit(1)
		# split info by sections
		self.split_sections()
		# info run read cycles
		self.infer_read_cycles()
		# extract samples by lane
		self.lane_info = self.load_section_table('BCLConvert_Data')
		# fix data type for column 'Lane'
		self.lane_info['Lane'] = self.lane_info['Lane'].astype('int32')
		# fix data type for column 'Sample_ID'
		self.lane_info['Sample_ID'] = self.lane_info['Sample_ID'].astype('string')
		# fill in OverrideCycles column in case it is not provided in sample sheet
		if 'OverrideCycles' not in self.lane_info.columns:
			self.lane_info['OverrideCycles'] = self.override_cycles
		# fill in index columns in case they are not provided in sample sheet
		if 'Index' not in self.lane_info.columns:
			self.lane_info['Index'] = ''
		if 'Index2' not in self.lane_info.columns:
			self.lane_info['Index2'] = ''
		# extract samples by project
		self.project_info = self.load_section_table('Cloud_Data')
		# fix data type for column 'Sample_ID'
		self.project_info['Sample_ID'] = self.project_info['Sample_ID'].astype('string')

	def format_override_cycles(self, tx):
		"""convert OverrideCycles to Description"""
		return '+'.join(['{}'.format(sum([int(k) for k in re.findall("[Y,I,U](\d+)", tc)])) for tc in tx.split(';')])

	def convert_v1_format(self):
		"""convert samplesheet back to v1 format"""
		logging.info('Convert samplesheet to v1 format.')
		# merge lane and project info tables
		self.v1_format = self.lane_info[['Lane','Sample_ID','Index','Index2','OverrideCycles']].merge(self.project_info[['Sample_ID','ProjectName']], on='Sample_ID', how='left')
		# reformat OverrideCycles
		self.v1_format['Description'] = self.v1_format['OverrideCycles'].apply(self.format_override_cycles)
		# rename columns
		self.v1_format.rename(columns={'Index':'index','Index2':'index2','ProjectName':'Sample_Project'}, inplace=True)
		# reorder lane and sample
		self.v1_format.sort_values(['Lane','Sample_Project'], inplace=True)
		# replace na index to empty str
		self.v1_format['index2'] = self.v1_format['index2'].replace({'na':''})
		# add missing columns
		self.v1_format['Sample_Name'] = ''
		self.v1_format['Sample_Plate'] = ''
		self.v1_format['Sample_Well'] = ''
		self.v1_format['Index_Plate_Well'] = ''
		self.v1_format['I7_Index_ID'] = ''
		self.v1_format['I5_Index_ID'] = ''

	def write_v1_format(self):
		"""write v1 format samplesheet to file"""
		if self.v1_format_file is None:
			return None
		logging.info('Write v1 format samplesheet to file.')
		# write headers
		with open(self.v1_format_file, 'w') as fout:
			fout.write("[Data],,,,,,,,,,,\n")
		self.v1_format.to_csv(self.v1_format_file, columns=self.v1_format_cols, sep=',', index=False, mode='a')

	def parse_fastq_files(self):
		"""parse fastq files and extract sample/lane/read info"""
		fastq_files = [x for x in os.listdir(self.demux_fastq_folder) if re.search("\.fastq\.gz$", x)]
		seq_info = [re.search("(.*?)_S(\d+)_L00(\d+)_([RI]+\d+)_001\.fastq\.gz", x).groups() for x in fastq_files]
		###print([x for x in seq_info if x[0] not in self.project_info['Sample_ID'].tolist()])
		return pd.DataFrame({'Sample_ID':[x[0] for x in seq_info], 'SN':[x[1] for x in seq_info], 'Lane':[x[2] for x in seq_info], \
			'Read':[x[3] for x in seq_info], 'Fastq_Folder':self.demux_fastq_folder, 'Fastq_File':fastq_files}).merge(\
			self.project_info[['Sample_ID','ProjectName']], on='Sample_ID', how='right')

	def parse_fastqc_report_files(self):
		"""parse fastqc report and extract sample info"""
		fastqc_report_folders = []
		fastqc_report_files = []
		fastqc_metrics_files = []
		samples = []
		for root, folder, files in os.walk(self.demux_bclconvert_folder):
			if re.search("fastqc$", root):
				# get sample id
				sid = os.path.basename(os.path.dirname(root))
				# report
				report_file = os.path.join(root, "report.html")
				if not os.path.exists(report_file):
					logging.error('Failed to locate fastqc report file: {}.'.format(report_file))
					#sys.exit(5)
					fastqc_report_files.append('-')
				else:
					fastqc_report_files.append("report.html")
				# metrics
				metrics_file = os.path.join(root, "{}.fastqc_metrics.csv".format(sid))
				if not os.path.exists(metrics_file):
					logging.error('Failed to locate fastqc metrics file: {}.'.format(metrics_file))
					#sys.exit(5)
					fastqc_metrics_files.append('-')
				else:
					fastqc_metrics_files.append("{}.fastqc_metrics.csv".format(sid))
				# folder
				fastqc_report_folders.append(root)
				# sample
				samples.append(sid)
		return pd.DataFrame({'Sample_ID':samples, 'Fastqc_Report_Folder':fastqc_report_folders, 'Fastqc_Report_File':fastqc_report_files, \
			'Fastqc_Metrics_File':fastqc_metrics_files, \
			'Fastqc_Report_File_rename':['fastqc.{}.html'.format(x) for x in samples]}).merge(\
			self.project_info[['Sample_ID','ProjectName']], on='Sample_ID', how='right')

	def scan_demux_files(self):
		"""scan demux folder and collect files needed for distribution"""
		if self.demux_folder is None:
			return None
		logging.info('Scan demux folder and collect files needed for distribution:')
		# demux stats file
		self.demux_stats_file = os.path.join(self.demux_folder, 'Data', 'Demux', 'Demultiplex_Stats.csv')
		# demux BCLConvert folder
		self.demux_bclconvert_folder = os.path.join(self.demux_folder, 'Data', 'BCLConvert')
		# demux fastq folder
		self.demux_fastq_folder = os.path.join(self.demux_folder, 'Data', 'BCLConvert', 'fastq')
		# check existence
		if not os.path.exists(self.demux_stats_file):
			logging.error('Failed to locate demux stats file: {}.'.format(self.demux_stats_file))
			sys.exit(4)
		if not os.path.exists(self.demux_bclconvert_folder):
			logging.error('Failed to locate demux BCLConvert folder: {}.'.format(self.demux_bclconvert_folder))
			sys.exit(5)
		if not os.path.exists(self.demux_fastq_folder):
			logging.error('Failed to locate demux fastq folder: {}.'.format(self.demux_fastq_folder))
			sys.exit(6)
		# extract demux stats table
		self.demux_stats = pd.read_table(self.demux_stats_file, sep=',', header=0, low_memory=False, dtype={'SampleID':'string'}).merge(\
			self.project_info[['Sample_ID','ProjectName']].rename(columns={'Sample_ID':'SampleID'}), on='SampleID', how='left')
		logging.info('Load demux stats file: {}.'.format(self.demux_stats.shape))
		# demux fastq files
		self.demux_fastq_files = self.parse_fastq_files()
		logging.info('Collect fastq files: {}.'.format(self.demux_fastq_files.shape))
		# demux fastqc report files
		self.demux_fastqc_report_files = self.parse_fastqc_report_files()
		logging.info('Collect fastqc report files: {}.'.format(self.demux_fastqc_report_files.shape))

	def initialize_output_folders(self):
		"""initialize output folder for each project"""
		if self.output_folder is None:
			return None
		if os.path.exists(self.output_folder):
			logging.info('Output folder already exists: {}.\nPlease delete the folder and rerun!'.format(self.output_folder))
			sys.exit(0)
		else:
			# create output folder
			try:
				os.mkdir(self.output_folder)
				# create folder for each project
				for project in self.project_info['ProjectName'].unique():
					os.mkdir(os.path.join(self.output_folder, project))
					os.mkdir(os.path.join(self.output_folder, project, 'Summary'))
					os.mkdir(os.path.join(self.output_folder, project, 'Fastqc'))
					os.mkdir(os.path.join(self.output_folder, project, 'MultiQC'))
			except OSError as error:
				logging.error(error)
				sys.exit(7)
			# script file
			self.arrange_fastq_script = os.path.join(self.output_folder, 'link_fastq.sh')
			self.arrange_fastqc_report_script = os.path.join(self.output_folder, 'link_fastqc_report.sh')
			self.run_multiqc_script = os.path.join(self.output_folder, 'run_multiqc.sh')
			self.run_multiqc_script_log = os.path.join(self.output_folder, 'run_multiqc.log')

	def subset_demux_stats(self):
		"""subset demux stats by project"""
		if self.output_folder is None:
			return None
		logging.info('Subset demux stats by project.')
		for project in self.project_info['ProjectName'].unique():
			# file name
			stats_filename = os.path.join(self.output_folder, project, 'Summary', 'Demultiplex_Stats.{}.html'.format(project))
			# write to file
			self.demux_stats[self.demux_stats['ProjectName'] == project].to_html(stats_filename, columns=self.demux_cols, \
				index=False, formatters={'# Reads': lambda x: '{:,}'.format(x)}, justify='left')

	def link_file_to_project(self, entry, folder_col, file_col, new_file_col, fs, target_folder, symbolic):
		"""link a file to a project"""
		foldername = entry[folder_col]
		filename = entry[file_col]
		filenewname = entry[new_file_col]
		# handle space in project name
		project = entry['ProjectName'].replace(' ', '\ ')
		outdir = os.path.join(self.output_folder, project, target_folder)
		if target_folder == '':
			outdir = os.path.join(self.output_folder, project)
		# make sure the file is available
		if filename != '-':
			fs.write('cd {}\n'.format(outdir))
			if symbolic:
				fs.write('ln -s {} {}\n'.format(os.path.join(foldername, filename), filenewname))
			else:
				fs.write('ln {} {}\n'.format(os.path.join(foldername, filename), filenewname))
			fs.write('\n')

	def arrange_fastq_files(self):
		"""arrange fastq files by project"""
		if self.output_folder is None:
			return None
		logging.info('Arrange fastq files by project.')
		with open(self.arrange_fastq_script, 'w') as fout:
			self.demux_fastq_files.apply(self.link_file_to_project, axis=1, args=('Fastq_Folder','Fastq_File','Fastq_File',fout,'',self.symbolic))
		# run script
		if not self.norun:
			try:
				subprocess.run(['sh', self.arrange_fastq_script], check=True)
			except subprocess.CalledProcessError as e:
				logging.error('Command {} failed with error {}.'.format(e.cmd, e.returncode))
				sys.exit(8)

	def arrange_fastqc_report_files(self):
		"""arrange fastqc report files by project"""
		if self.output_folder is None:
			return None
		logging.info('Arrange fastqc report files by project.')
		with open(self.arrange_fastqc_report_script, 'w') as fout:
			self.demux_fastqc_report_files.apply(self.link_file_to_project, axis=1, args=('Fastqc_Report_Folder','Fastqc_Report_File','Fastqc_Report_File_rename',fout,'Fastqc',self.symbolic))
		# run script
		if not self.norun:
			try:
				subprocess.run(['sh', self.arrange_fastqc_report_script], check=True)
			except subprocess.CalledProcessError as e:
				logging.error('Command {} failed with error {}.'.format(e.cmd, e.returncode))
				sys.exit(9)

	def run_multiqc(self):
		"""run multiqc to combine fastqc results from multiple samples by project"""
		if self.output_folder is None:
			return None
		logging.info('Running multiqc to combine fastqc results by project.')
		if self.skip_multqc is not None:
			logging.info('Skip multiqc for projects: {}'.format(self.skip_multqc))
			self.skip_multqc = self.skip_multqc.split(';')
		else:
			self.skip_multqc = []
		with open(self.run_multiqc_script, 'w') as fs:
			project_list = self.demux_fastqc_report_files['ProjectName'].unique()
			for project in project_list:
				if project in self.skip_multqc:
					continue
				# subset by project
				td = self.demux_fastqc_report_files[self.demux_fastqc_report_files['ProjectName'] == project]
				# prepare codes
				fs.write('echo {}\n'.format(project))
				fs.write('{} \\\n'.format(self.multiqc))
				# handle space in project name
				project = project.replace(' ', '\ ')
				fs.write(''.join(['\t{} \\\n'.format(os.path.join(folder, filename)) for folder, filename in zip(td['Fastqc_Report_Folder'].tolist(), td['Fastqc_Metrics_File'].tolist())]))
				fs.write('\t-o {} -n {} \n\n'.format(os.path.join(self.output_folder, project, 'MultiQC'), 'multiqc_report.{}.html'.format(project)))
		# run script
		if not self.norun:
			try:
				result = subprocess.run(['sh', self.run_multiqc_script], check=True, capture_output=True)
				with open(self.run_multiqc_script_log, 'w') as flog:
					flog.write(result.stdout.decode())
					flog.write(result.stderr.decode())
			except subprocess.CalledProcessError as e:
				logging.error('Command {} failed with error {}.'.format(e.cmd, e.returncode))
				sys.exit(10)

	def generate_demux_summary(self):
		"""generate an overall demux summary by project"""
		if self.output_folder is None:
			return None
		logging.info('Generate overall demux summary by project.')
		# count reads for each project
		self.demux_summary = self.demux_stats.groupby('ProjectName')['# Reads'].sum().to_frame(name='# Total Reads').reset_index()
		# output file name
		summary_filename = os.path.join(self.output_folder, self.demux_summary_file)
		# write to file
		self.demux_summary.to_html(summary_filename, index=False, formatters={'# Total Reads': lambda x: '{:,}'.format(x)}, justify='left')

# function
def get_arguments():
	"""fetch command line arguments."""
	parser = ArgumentParser(description="""Collect demux data processed by Illumina BCL Convert, pack fastq files and prepare summary reports for distribution.""",
                            prog='arrange_fastq.py')
	parser.add_argument("-v", "--version", action="version", version='%(prog)s v0.3')
	parser.add_argument("-i", "--samplesheet", nargs="?", required=True, help="samplesheet made on basespace", metavar="samplesheet_file", dest="samplesheetfile")
	parser.add_argument("-l", "--log", nargs="?", default="arrange_fastq.log", help="log file", metavar="log_file", dest="logfile")
	parser.add_argument("-c", "--v1format", nargs="?", help="convert samplesheet to v1 format", dest="v1formatfile")
	parser.add_argument("-d", "--demuxfolder", nargs="?", help="demux folder", dest="demuxfolder")
	parser.add_argument("-o", "--outfolder", nargs="?", help="output folder containing re-arranged demux files", dest="outfolder")
	parser.add_argument("-s", "--symbolic", action='store_true', default=False, help="make symbolic links instead of hard links", dest="symbolic")
	parser.add_argument("-n", "--no-run-script", action='store_true', default=False, help="generate link-file shell scripts without running them", dest="norun")
	parser.add_argument("-r", "--run-summary", nargs="?", default="demux_summary.html", help="an overall demux summary by project", metavar="demux_summary_file", dest="sumfile")
	parser.add_argument("-q", "--multiqc", nargs="?", default="multiqc", help="multiqc excutable", dest="multiqc")
	parser.add_argument("-k", "--skip-multiqc", nargs="?", help="list projects that you want to skip multiqc, separate projects by comma(;)", dest="skipmultqc")
	return parser.parse_args()

def setup_logging(logfile, level=logging.INFO):
	"""set up logging"""
	# prepare loggings
	log_formatter = logging.Formatter('%(levelname)s: %(message)s')
	root_logger = logging.getLogger()
	root_logger.setLevel(level)

	# logging to file
	file_handler = logging.FileHandler(logfile)
	file_handler.setFormatter(log_formatter)
	root_logger.addHandler(file_handler)

	# logging to stdout
	console_handler = logging.StreamHandler(sys.stdout)
	console_handler.setFormatter(log_formatter)
	root_logger.addHandler(console_handler)

	# another common way to set up logging, but of less flexibility
	#logging.basicConfig(level=logging.DEBUG, format='%(levelname)s: %(message)s')
	#logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
	return root_logger

def main():
	"""call me to get started!"""
	# load arguments
	args = get_arguments()
	# set up logging
	root_logger = setup_logging(args.logfile, level=logging.INFO)
	#root_logger = setup_logging(args.logfile, level=logging.DEBUG)
	#root_logger = setup_logging(args.logfile, level=logging.INFO)
	# create a MyBCLConvert object
	d = MyBCLConvert(args)
	# extract information from sample sheet file
	d.get_sample_info()
	##print(d.sections)
	##print(d.lane_info.shape)
	##print(d.lane_info.head())
	##print(d.lane_info.tail())
	##print(d.project_info.shape)
	##print(d.project_info.head())
	##print(d.project_info.tail())
	# convert samplesheet to v1 format
	##print(d.v1_format_file)
	d.convert_v1_format()
	d.write_v1_format()
	# scan demux files
	d.scan_demux_files()
	##print(d.demux_stats.head())
	##print(d.demux_fastq_files.head())
	##print(d.demux_fastqc_report_files.head())
	##print(d.project_info['ProjectName'].unique())
	##print(d.demux_summary.head())
	# re-arrange demux files
	d.initialize_output_folders()
	d.subset_demux_stats()
	d.arrange_fastq_files()
	d.arrange_fastqc_report_files()
	d.run_multiqc()
	# generate overall demux summary by project
	d.generate_demux_summary()

# main
if __name__ == '__main__':
	main()
