#!/usr/bin/env python
# collect_demux_info.py
# collect demux info from sequencing runs, for records
#

import sys
import gzip
import logging

import pandas as pd
import xml.etree.ElementTree as ET
import os.path

from os import getcwd
from os import listdir
from os import walk
from re import findall
from re import search
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter

# classes
# MyDemuxUnit: a single demux unit:
#     a set of samples saved in a separate demux folder, 
#     which usually has the same iLab request id & sequencing type/length
class MyDemuxUnit:
    """Class for saving info of a single demux unit"""

    def __init__(self, fastq_folder, flowcell_id, report_file=None, ilab=-1, project=None, seqtype=None, read_1_len=-1, read_2_len=-1):
        """class constructor"""
        # fastq data folder
        self.fastq_folder = fastq_folder
        # flowcell id
        self.fcid = flowcell_id
        # report file
        self.report_file = report_file
        # iLab request
        self.ilab = ilab
        # project name
        self.project = project
        # sample demux summary
        self.summary = None
        # sequencing type (SR, PE)
        self.seqtype = seqtype
        # read length
        self.read_1_len = read_1_len
        self.read_2_len = read_2_len
        # columns to include
        self.columns = ['Sample','Barcode sequence', 'PF Clusters', 'Yield (Mbases)', '% >= Q30bases']
        # a valid demux unit for output?
        self.valid = True

    def infer_ilab(self):
        """guest iLab id and project name"""
        # eg: Loda-MJ-10557_2021_06_15
        #     Diaz-Meco-MADM-10830_2021_07_19
        #     Delia
        #     Delia_2021_10_15
        if self.ilab == -1:
            folder_name = os.path.basename(self.fastq_folder)
            tmatches_A = findall('(\d+)', folder_name)
            tmatches_B = search('(.*?)_\d+_\d+_\d+', folder_name)
            if tmatches_A:
                self.ilab = int(tmatches_A[0])
                self.project = folder_name[:folder_name.find(tmatches_A[0])-1]
                logging.debug('Infer iLab id: {}'.format(self.ilab))
                logging.debug('Infer project info: {}'.format(self.project))
            else:
                # some project does not have a iLab id, so leave it as '-1' while still keep this project
                ###self.valid = False
                # get project name
                self.project = folder_name
                if tmatches_B:
                    self.project = tmatches_B.groups()[0]
                logging.warning('Failed to infer iLab id: {}'.format(self.fastq_folder))
                logging.debug('Infer project info: {}'.format(self.project))

    def infer_report_file(self):
        """guess demux report html file if not provided"""
        # report file already exists?
        if self.report_file is None:
            # the 'Reports' folder should be in the same level as the fastq folder and 
            # they share the same parent folder, not necessarily named 'Unaligned' (e.g. 10X runs)
            tbasefolder = os.path.dirname(self.fastq_folder)
            # get project folder
            # e.g. Nam-SR-11103_2021_09_01_part2 --> Nam-SR-11103
            #      Nam-SR-11103 (before data distribution, no need to change)
            #      Delia (no need to change)
            tpjtfolder = os.path.basename(self.fastq_folder)
            tpat = search('^(.*?)_\d+_\d+_\d+', tpjtfolder)
            if tpat:
                tpjtfolder = tpat.groups()[0]
            # eg. Nam-SR-11103_2021_09_01_part2
            # guess report html file
            tsumfile = os.path.join(tbasefolder,'Reports','html',self.fcid,tpjtfolder,'all','all','laneBarcode.html')
            # check file existence
            if os.path.exists(tsumfile):
                self.report_file = tsumfile
                logging.debug('Infer demux report html file: {}'.format(self.report_file))
            else:
                self.valid = False
                logging.warning('Failed to locate demux report html file: {}'.format(tsumfile))

    def guess_read_length(self, fastq_file, num_reads=100):
        """guess read length for a given fastq file"""
        # sample the first 'num_reads' reads
        k = 0
        rlen = 0
        with gzip.open(fastq_file, 'r') as fin:
            # sample the first 'num_reads' reads
            for line in fin:
                if k % 4 == 1:# read sequence
                    clen = len(line.strip())
                    if clen > rlen:
                        rlen = clen
                k += 1
                if k > 4 * num_reads:
                    break
        return rlen

    def infer_seq_type(self):
        """guess sequencing type and read length based on fastq files"""
        if self.seqtype is None or self.read_1_len is None:
            # fetch available fastq file
            tfqs = []
            for root, directories, files in walk(self.fastq_folder):
                tfqs.extend([os.path.join(root, x) for x in files if search('fastq\.gz$',x) and not search('^Undetermined_',x)])
            # fastq files found?
            if len(tfqs) > 0:
                # read 1/2 fastq
                tr1s = [x for x in tfqs if search('R1_001\.fastq\.gz', x)]
                tr2s = [x for x in tfqs if search('R2_001\.fastq\.gz', x)]
                # SR or PE
                if len(tr2s) == 0:
                    self.seqtype = 'SR'
                else:
                    self.seqtype = 'PE'
                # read length
                self.read_1_len = self.guess_read_length(tr1s[0])
                if len(tr2s) > 0:
                    self.read_2_len = self.guess_read_length(tr2s[0])
                logging.debug('Infer sequencing type: {}'.format(self.seqtype))
                logging.debug('Infer read 1 length: {}'.format(self.read_1_len))
                if len(tr2s) > 0:
                    logging.debug('Infer read 2 length: {}'.format(self.read_2_len))

    def load_report(self):
        """extract demux summary from report file"""
        if self.report_file is not None:
            sumtables = pd.read_html(self.report_file, match='Sample')
            # double check on detected tables
            if len(sumtables) != 1:
                logging.warning("More than one table with 'Sample' column detected in demux summary file {}. Use the first table.".format(self.report_file))
            else:
                self.summary = sumtables[0]
            # double check if all required columns exist
            for col in self.columns:
                if col not in self.summary.columns:
                    self.valid = False
                    logging.warning("Column '{}' cannot be found in the demux summary file {}.".format(col, self.report_file))
            # replace 'Sample_' in the sample name
            self.summary['Sample'] = self.summary['Sample'].astype('str')
            self.summary['Sample'] = self.summary['Sample'].str.replace(r'^Sample_','',regex=True)
            # remove any Sample == 'Undetermined' & Barcode sequence == 'unknown'
            filt1 = self.summary['Sample'] == 'Undetermined'
            filt2 = self.summary['Barcode sequence'] != 'unknown'
            self.summary = self.summary[~ (filt1 & filt2)]

    def infer_all(self):
        """infer all available information"""
        self.infer_ilab()
        self.infer_report_file()
        self.infer_seq_type()
        self.load_report()

    def prepare_table(self):
        """combine information into a dataframe table"""
        # a valid demux unit?
        if not self.valid:
            logging.warning('Not a valid demux project: {}'.format(self.fastq_folder))
            return None
        # get a copy of summary for output purpose
        mytable = self.summary[self.columns].copy()
        # add other information
        mytable.insert(0, 'iLab', [self.ilab] * mytable.shape[0])
        mytable.insert(1, 'project', [self.project] * mytable.shape[0])
        mytable.insert(2, 'seqtype', [self.seqtype] * mytable.shape[0])
        mytable.insert(3, 'read1len', [self.read_1_len] * mytable.shape[0])
        mytable.insert(4, 'read2len', [self.read_2_len] * mytable.shape[0])
        return mytable

    def to_file(self, outfile):
        """combine information and write to file"""
        # a valid demux unit?
        if not self.valid:
            logging.warning('Not a valid demux project: {}'.format(self.fastq_folder))
            return None
        # combine information into a table
        mytable = self.prepare_table()
        # write to file
        if search('\.xlsx$', outfile):
            mytable.to_excel(outfile, index=False)
        elif search('\.csv$', outfile):
            mytable.to_csv(outfile, sep=',', index=False)
        elif search('\.tsv$|\.txt$', outfile):
            mytable.to_csv(outfile, sep='\t', index=False)
        else:
            logging.warning('Unsupported output file extension: {}'.format(outfile.split('.')[-1]))
            return None
        logging.info('write MyDemuxUnit to file: {}'.format(outfile))

    def __repr__(self):
        """print info for a demux unit (project)"""
        to_print = []
        to_print.append('fastq folder: {}'.format(self.fastq_folder))
        to_print.append('flowcell id: {}'.format(self.fcid))
        to_print.append('report file: {}'.format(self.report_file))
        to_print.append('iLab id: {}'.format(self.ilab))
        to_print.append('project name: {}'.format(self.project))
        to_print.append('sequence type: {}'.format(self.seqtype))
        to_print.append('read 1 length: {}'.format(self.read_1_len))
        to_print.append('read 2 length: {}'.format(self.read_2_len))
        to_print.append('valid?: {}'.format(self.valid))
        return '\n'.join(to_print)

# MyDemuxRun: a sequencing run
class MyDemuxRun:
    """Class for saving info of an entire demux run"""

    def __init__(self, run_folder, date=None, flowcell_id=None, platform=None, seqtype=None, read_1_len=-1, read_2_len=-1, index_1_len=-1, index_2_len=-1):
        """class constructor"""
        # run folder
        self.run_folder = run_folder
        # sequencing date
        self.date = date
        # flowcell id
        self.fcid = flowcell_id
        # sequencing platform
        self.platform = platform
        # run sequencing type
        self.seqtype = seqtype
        # run read length
        self.read_1_len = read_1_len
        self.read_2_len = read_2_len
        self.index_1_len = index_1_len
        self.index_2_len = index_2_len
        # projects in this run
        self.units = []
        # a valid demux run for output?
        self.valid = True

    def infer_platform(self):
        """infer platform from run folder"""
        if self.platform is None:
            if 'hiseq' in self.run_folder.lower():
                self.platform = 'hiseq'
            elif 'nextseq500' in self.run_folder.lower():
                self.platform = 'nextseq500'
            elif 'nextseq2000' in self.run_folder.lower():
                self.platform = 'nextseq2000'
            elif 'novaseq' in self.run_folder.lower():
                self.platform = 'novaseq'
            else:
                self.valid = False
                logging.warning('Failed to infer platform for {}'.format(self.platform))
            logging.debug('Infer platform: {}'.format(self.platform))

    def verify_seq_date(self, date):
        """a date (e.g. 201208) is valid?"""
        # include non-numeric characters?
        if search('\D+', date):
            return False
        # month
        if date[2] not in ['0','1']:
            return False
        # day
        if date[4] not in ['0','1','2','3']:
            return False
        return True

    def infer_seq_date(self):
        """guess date of sequencing"""
        if self.date is None:
            tdate = os.path.basename(self.run_folder).split('_')[0]
            # double check on date, make sure it is a valid date
            if not self.verify_seq_date(tdate):
                self.valid = False
                logging.warning('Failed to infer sequencing date for {}'.format(self.run_folder))
            else:
                self.date = '20{}-{}-{}'.format(tdate[:2],tdate[2:4],tdate[4:])
                logging.debug('Infer sequencing date: {}'.format(self.date))

    def infer_flowcell_id(self):
        """guess flowcell id"""
        if self.fcid is None:
            self.fcid = os.path.basename(self.run_folder).split('_')[-1]
            if self.platform not in ['nextseq2000','miseq']:
                self.fcid = self.fcid[1:]# the first letter refers to flowcell A or B
            logging.debug('Infer flowcell id: {}'.format(self.fcid))

    def infer_read_length(self):
        """extract sequenced read length from runparameter.xml"""
        # already known, do nothing
        if self.read_1_len > 0:
            return None
        # default setting for novaseq
        # run parameter file
        run_para_file = '{}/RunParameters.xml'.format(self.run_folder)
        if self.platform == 'hiseq':
            run_para_file = '{}/runParameters.xml'.format(self.run_folder)
        # read length labels
        r1id = 'Read1NumberOfCycles'
        r2id = 'Read2NumberOfCycles'
        i1id = 'IndexRead1NumberOfCycles'
        i2id = 'IndexRead2NumberOfCycles'
        if self.platform == 'hiseq':
            r1id = 'Read1'
            r2id = 'Read2'
            i1id = 'IndexRead1'
            i2id = 'IndexRead2'
        elif self.platform == 'nextseq500':
            r1id = 'Read1'
            r2id = 'Read2'
            i1id = 'Index1Read'
            i2id = 'Index2Read'
        elif self.platform == 'nextseq2000':
            r1id = 'Read1'
            r2id = 'Read2'
            i1id = 'Index1'
            i2id = 'Index2'
        # run parameter file exists?
        if not os.path.exists(run_para_file):
            self.valid = False
            logging.warning('Unable to read the run parameter file: {}'.format(run_para_file))
            return None
        # parse xml file
        tree = ET.parse(run_para_file)
        # get root
        root = tree.getroot()
        # extract read length info
        for x in root.iter(r1id):
            self.read_1_len = int(x.text)
        for x in root.iter(r2id):
            self.read_2_len = int(x.text)
        for x in root.iter(i1id):
            self.index_1_len = int(x.text)
        for x in root.iter(i2id):
            self.index_2_len = int(x.text)
        # do we have length for all reads?
        if self.read_1_len == -1:
            self.valid = False
            logging.warning('Failed to get read 1 length.')
        elif self.read_2_len == -1:
            self.valid = False
            logging.warning('Failed to get read 2 length.')
        elif self.index_1_len == -1:
            self.valid = False
            logging.warning('Failed to get index 1 length.')
        elif self.index_2_len == -1:
            self.valid = False
            logging.warning('Failed to get index 2 length.')
        else:
            logging.debug('Infer read 1 length: {}'.format(self.read_1_len))
            logging.debug('Infer read 2 length: {}'.format(self.read_2_len))
            logging.debug('Infer index 1 length: {}'.format(self.index_1_len))
            logging.debug('Infer index 2 length: {}'.format(self.index_2_len))

    def infer_seq_type(self):
        """infer sequencing type based on read length"""
        if self.seqtype is None:
            if self.read_2_len > 0:
                self.seqtype = 'PE'
            else:
                self.seqtype = 'SR'
            logging.debug('Infer sequencing type: {}'.format(self.seqtype))

    def extract_demux_units(self):
        """extract demux information for each project"""
        # screen for fastq folders under run folder
        # assume fastq stored under folder named 'Unaligned', but not necessarily directly under it (e.g. 10X runs)
        candidate_folders = [x for x in listdir(self.run_folder) if search('Unaligned',x)]
        # scan for folders containing fastq files (demux unit)
        temp_folders = []
        for folder in candidate_folders:
            for root, directories, files in walk(os.path.join(self.run_folder, folder)):
                for tfile in files:
                    if search('\.fastq\.gz$', tfile) and not search('^Undetermined_', tfile):
                        # there are cases where fastq files were saved in a separate folder per sample
                        if not search('.*?\-.*?\-\d+', os.path.basename(root)):
                            temp_folders.append(os.path.dirname(root))
                        else:
                            temp_folders.append(root)
                        break
        fastq_folders = list(set(temp_folders))
        #logging.debug('\n'.join(fastq_folders))
        # process each demux unit and store it as a MyDemuxUnit object
        for folder in fastq_folders:
            logging.debug(folder)
            tunit = MyDemuxUnit(folder, self.fcid)
            tunit.infer_all()
            if tunit.valid:
                self.units.append(tunit)
        # at least one valid unit collected?
        if not self.units:
            self.valid = False
            logging.warning('Not a valid demux run since no valid units (projects) are found: {}'.format(self.run_folder))

    def infer_all(self):
        """infer all available information"""
        self.infer_platform()
        self.infer_seq_date()
        self.infer_flowcell_id()
        self.infer_read_length()
        self.infer_seq_type()
        self.extract_demux_units()

    def prepare_table(self):
        """combine information into a dataframe table"""
        # a valid demux run?
        if not self.valid:
            return None
        # merge information from each demux unit
        mytable = pd.concat([x.prepare_table() for x in self.units])
        # add additional information
        # platform
        mytable.insert(2, 'platform', [self.platform] * mytable.shape[0])
        # flowcell id
        mytable.insert(3, 'flowcell_id', [self.fcid] * mytable.shape[0])
        # sequencing date
        mytable.insert(4, 'date', [self.date] * mytable.shape[0])
        mytable['date'] = pd.to_datetime(mytable['date'], format="%Y-%m-%d")
        # run sequencing type
        mytable.insert(mytable.shape[1], 'run_seqtype', [self.seqtype] * mytable.shape[0])
        # run read length
        mytable.insert(mytable.shape[1], 'run_read1len', [self.read_1_len] * mytable.shape[0])
        mytable.insert(mytable.shape[1], 'run_read2len', [self.read_2_len] * mytable.shape[0])
        mytable.insert(mytable.shape[1], 'run_index1len', [self.index_1_len] * mytable.shape[0])
        mytable.insert(mytable.shape[1], 'run_index2len', [self.index_2_len] * mytable.shape[0])
        #logging.debug(mytable.shape)
        #logging.debug(mytable.head(2))
        return mytable

    def to_file(self, outfile):
        """combine information and write to file"""
        # a valid demux run?
        if not self.valid:
            return None
        # combine information into a table
        mytable = self.prepare_table()
        # write to file
        if search('\.xlsx$', outfile):
            mytable.to_excel(outfile, index=False)
        elif search('\.csv$', outfile):
            mytable.to_csv(outfile, sep=',', index=False)
        elif search('\.tsv$|\.txt$', outfile):
            mytable.to_csv(outfile, sep='\t', index=False)
        else:
            logging.warning('Unsupported output file extension: {}'.format(outfile.split('.')[-1]))
            return None
        logging.info('write MyDemuxRun to file: {}'.format(outfile))

    def __repr__(self):
        """print info for a demux run"""
        to_print = []
        to_print.append('run folder: {}'.format(self.run_folder))
        to_print.append('platform: {}'.format(self.platform))
        to_print.append('flowcell id: {}'.format(self.fcid))
        to_print.append('date: {}'.format(self.date))
        to_print.append('sequence type: {}'.format(self.seqtype))
        to_print.append('read 1 length: {}'.format(self.read_1_len))
        to_print.append('read 2 length: {}'.format(self.read_2_len))
        to_print.append('index 1 length: {}'.format(self.index_1_len))
        to_print.append('index 2 length: {}'.format(self.index_2_len))
        to_print.append('# demux units: {}'.format(len(self.units)))
        to_print.append('valid?: {}'.format(self.valid))
        return '\n'.join(to_print)

# MyDemuxFolder: a folder under which multiple sequencing runs lie
class MyDemuxFolder:
    """Class for saving info of runs under a given folder"""

    def __init__(self, server_folder):
        """class constructor"""
        # server folder
        self.server_folder = server_folder
        # sequencing runs inside this folder
        self.runs = []
        # a valid demux folder for output?
        self.valid = True

    def extract_demux_runs(self):
        """extract demux information for all available runs inside the folder"""
        logging.info('[Folder]\t{}'.format(self.server_folder))
        # screen for all available sequencing runs
        # e.g. 210901_A00814_0481_AHJ5T2DRXY
        run_folders = [x for x in listdir(self.server_folder) if search('\d+_[A-Z\d]+_\d+_[A-Za-z\d]+',x)]
        # process each demux run folder and store it as a MyDemuxRun object
        for folder in run_folders:
            logging.info('[Run]\t{}'.format(folder))
            trun = MyDemuxRun(os.path.join(self.server_folder,folder))
            trun.infer_all()
            if trun.valid:
                self.runs.append(trun)
            #logging.debug(trun)
        # at least one valid demux run?
        if not self.runs:
            self.valid = False
            logging.warning('Not a valid demux folder since no valid runs are found: {}'.format(self.server_folder))

    def prepare_table(self):
        """combine information into a dataframe table"""
        # a valid demux folder?
        if not self.valid:
            return None
        # merge information from each demux run
        mytable = pd.concat([x.prepare_table() for x in self.runs])
        #logging.debug(mytable.shape)
        #logging.debug(mytable.head(2))
        return mytable

    def to_file(self, outfile):
        """combine information and write to file"""
        # a valid demux folder?
        if not self.valid:
            return None
        # combine information into a table
        mytable = self.prepare_table()
        # write to file
        if search('\.xlsx$', outfile):
            mytable.to_excel(outfile, index=False)
        elif search('\.csv$', outfile):
            mytable.to_csv(outfile, sep=',', index=False)
        elif search('\.tsv$|\.txt$', outfile):
            mytable.to_csv(outfile, sep='\t', index=False)
        else:
            logging.warning('Unsupported output file extension: {}'.format(outfile.split('.')[-1]))
            return None
        logging.info('write MyDemuxFolder to file: {}'.format(outfile))

    def __repr__(self):
        """print info for a demux folder"""
        to_print = []
        to_print.append('server folder: {}'.format(self.server_folder))
        to_print.append('# demux runs: {}'.format(len(self.runs)))
        to_print.append('valid?: {}'.format(self.valid))
        return '\n'.join(to_print)

# MyDemuxAuto: automatically screen multiple folders and collect demux info
class MyDemuxAuto:
    """Class for saving info on multiple server locations"""

    def __init__(self, server_folder_list):
        """class constructor"""
        # a list of server folders to screen
        self.server_folder_list = server_folder_list
        # demux info for all available server folders
        self.folders = []
        # a valid demux auto for output?
        self.valid = True

    def extract_demux_folders(self):
        """extract demux information for all available folders in the list"""
        for folder in self.server_folder_list:
            tfolder = MyDemuxFolder(folder)
            tfolder.extract_demux_runs()
            if tfolder.valid:
                self.folders.append(tfolder)
            #logging.debug(tfolder)
        # at least one valid demux folder?
        if not self.folders:
            self.valid = False
            logging.warning('Not a valid demux auto since no valid folders are found: \n{}'.format('\n'.join(self.server_folder_list)))

    def prepare_table(self):
        """combine information into a dataframe table"""
        # a valid demux auto?
        if not self.valid:
            return None
        # merge information from all demux folders
        mytable = pd.concat([x.prepare_table() for x in self.folders])
        #logging.debug(mytable.shape)
        #logging.debug(mytable.head(2))
        return mytable

    def to_file(self, outfile):
        """combine information and write to file"""
        # a valid demux auto?
        if not self.valid:
            return None
        # combine information into a table
        mytable = self.prepare_table()
        # write to file
        if search('\.xlsx$', outfile):
            mytable.to_excel(outfile, index=False)
        elif search('\.csv$', outfile):
            mytable.to_csv(outfile, sep=',', index=False)
        elif search('\.tsv$|\.txt$', outfile):
            mytable.to_csv(outfile, sep='\t', index=False)
        else:
            logging.warning('Unsupported output file extension: {}'.format(outfile.split('.')[-1]))
            return None
        logging.info('write MyDemuxAuto to file: {}'.format(outfile))

    def __repr__(self):
        """print info for a demux auto"""
        to_print = []
        to_print.append('# demux folders: {}'.format(len(self.folders)))
        to_print.append('valid?: {}'.format(self.valid))
        return '\n'.join(to_print)

# functions
def test_MyDemuxUnit():
    #du = MyDemuxUnit('/gc7-data/NovaSeq6000/211015_A00814_0503_AHL5G2DSX2/Unaligned_1/Guzman-NMT-11319_2021_10_15','HL5G2DSX2')
    du = MyDemuxUnit('/gc7-data/NovaSeq6000/211015_A00814_0503_AHL5G2DSX2/Unaligned_2/Mason-ND-11312_2021_10_15','HL5G2DSX2')
    #du.infer_ilab()
    #du.infer_report_file()
    #du.infer_seq_type()
    #du.load_report()
    du.infer_all()
    du.to_file('/tmp/test.txt')
    print(du)
    #print(du.summary.columns)

def test_MyDemuxRun():
    dr = MyDemuxRun('/gc7-data/NovaSeq6000/211015_A00814_0503_AHL5G2DSX2')
    #dr = MyDemuxRun('/scratch/seq_data/NovaSeq6000/211014_A00814_0502_BHL5H3DSX2')
    #dr.infer_platform()
    #dr.infer_seq_date()
    #dr.infer_flowcell_id()
    #dr.infer_read_length()
    #dr.infer_seq_type()
    #dr.extract_demux_units()
    #dr.prepare_table()
    dr.infer_all()
    dr.to_file('/tmp/test.run.txt')
    print(dr)
    #print(len(dr.units))
    #print(dr.units[0])

def test_MyDemuxFolder():
    df = MyDemuxFolder('/scratch/seq_data/NovaSeq6000')
    #print(df.server_folder)
    df.extract_demux_runs()
    df.prepare_table()
    df.to_file('/tmp/test.folder.txt')
    print(df)

def test_MyDemuxAuto():
    da = MyDemuxAuto(['/scratch/seq_data/NovaSeq6000','/data/seq/NovaSeq6000'])
    #print(da.server_folder_list)
    da.extract_demux_folders()
    da.prepare_table()
    da.to_file('/tmp/test.auto.txt')

def main():
    #logging.basicConfig(level=logging.DEBUG, format='%(levelname)s: %(message)s')
    #logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    
    # prepare loggings
    log_formatter = logging.Formatter('%(levelname)s: %(message)s')
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.INFO)

    # logging to file
    file_handler = logging.FileHandler('/tmp/test.auto.log')
    file_handler.setFormatter(log_formatter)
    root_logger.addHandler(file_handler)

    # logging to stdout
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(log_formatter)
    root_logger.addHandler(console_handler)

    #test_MyDemuxUnit()
    #test_MyDemuxRun()
    #test_MyDemuxFolder()
    test_MyDemuxAuto()

# main
if __name__ == '__main__':
    main()

