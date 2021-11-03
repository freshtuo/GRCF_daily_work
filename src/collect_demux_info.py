#!/usr/bin/env python
# collect_demux_info.py
# collect demux info from sequencing runs, for records
#

import sys
import os
import gzip
import logging

import pandas as pd
import xml.etree.ElementTree as ET

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

    def __init__(self, fastq_folder, flowcell_id, report_file=None, ilab=-1, project=None, seqtype=None, read_1_len=-1, read_2_len=-1, detail_table=None):
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
        # detailed data table
        self.detail = detail_table
        # columns to include
        self.columns = ['Sample','Barcode sequence','num_lanes','PF Clusters','Yield (Mbases)','% >= Q30bases']
        # a valid demux unit for output?
        self.valid = True

    def infer_ilab(self):
        """guest iLab id and project name"""
        # assume 4-5 digit iLab id
        # eg: Loda-MJ-10557_2021_06_15         --> Loda-MJ & 10557
        #     Diaz-Meco-MADM-10830_2021_07_19  --> Diaz-Meco-MADM & 10830
        #     CS-10045_2021_03_25              --> CS & 10045
        #     Rebecca5830_2018_04_05           --> Rebecca & 5830
        #     Ning6652_rerun_2019_03_13
        #     Kaifang6912
        #     Delia                            --> Delia & -1
        #     Delia_2021_10_15                 --> Delia & -1
        #     Mason-DB-8850-2                  --> Mason-DB & 8850
        #     Nuria_2019_02_08_lane1           --> Nuria & -1
        #     Shawon_dropseq_062216            --> 
        #     Shawon_Suture_S4                 --> 
        #     Asaf07062016
        if self.ilab == -1:
            folder_name = os.path.basename(self.fastq_folder)
            tmatches_A = search('(\D+)\-([\d]{4,5})_[\d]{4}_[\d]{2}_[\d]{2}', folder_name)# (XXXX)-(10830)_2021_07_19
            tmatches_B = search('(\D+)([\d]{4,5})_[\d]{4}_[\d]{2}_[\d]{2}', folder_name)# (XXXX)(5830)_2018_04_05
            tmatches_C = search('(\D+)([\d]{4,5})_.*?_[\d]{4}_[\d]{2}_[\d]{2}', folder_name)# (Ning)(6652)_rerun_2019_03_13
            tmatches_D = search('(\D+)-([\d]{4,5})$', folder_name)# (Mason-DB)-(8850)
            tmatches_E = search('^(\D+)-([\d]{4,5})\D+', folder_name)# (Mason-DB)-(8850)-2
            tmatches_F = search('(\D+)([\d]{4,5})$', folder_name)# (Kaifang)(6912)
            tmatches_G = search('(\D+)_([\d]{4,5})$', folder_name)# (XXXX)_(5002)
            tmatches_H = search('(\D+)_[\d]{4}_[\d]{2}_[\d]{2}', folder_name)# (XXXX)_2021_10_15 or Nuria_2019_02_08_lane1
            tmatches_I = search('([A-Za-z]+)[_\-]*[\d]{6,8}$', folder_name)# (XXXX)_211015 or Nuria_20190208
            if tmatches_A:
                self.project, self.ilab = tmatches_A.groups()
                logging.debug('MyDemuxUnit: Infer iLab id: {}'.format(self.ilab))
                logging.debug('MyDemuxUnit: Infer project info: {}'.format(self.project))
            elif tmatches_B:
                self.project, self.ilab = tmatches_B.groups()
                logging.debug('MyDemuxUnit: Infer iLab id: {}'.format(self.ilab))
                logging.debug('MyDemuxUnit: Infer project info: {}'.format(self.project))
            elif tmatches_C:
                self.project, self.ilab = tmatches_C.groups()
                logging.debug('MyDemuxUnit: Infer iLab id: {}'.format(self.ilab))
                logging.debug('MyDemuxUnit: Infer project info: {}'.format(self.project))
            elif tmatches_D:
                self.project, self.ilab = tmatches_D.groups()
                logging.debug('MyDemuxUnit: Infer iLab id: {}'.format(self.ilab))
                logging.debug('MyDemuxUnit: Infer project info: {}'.format(self.project))
            elif tmatches_E:
                self.project, self.ilab = tmatches_E.groups()
                logging.debug('MyDemuxUnit: Infer iLab id: {}'.format(self.ilab))
                logging.debug('MyDemuxUnit: Infer project info: {}'.format(self.project))
            elif tmatches_F:
                self.project, self.ilab = tmatches_F.groups()
                logging.debug('MyDemuxUnit: Infer iLab id: {}'.format(self.ilab))
                logging.debug('MyDemuxUnit: Infer project info: {}'.format(self.project))
            elif tmatches_G:
                self.project, self.ilab = tmatches_G.groups()
                logging.debug('MyDemuxUnit: Infer iLab id: {}'.format(self.ilab))
                logging.debug('MyDemuxUnit: Infer project info: {}'.format(self.project))
            elif tmatches_H:
                self.project = tmatches_H.groups()[0]
                logging.warning('MyDemuxUnit: Failed to infer iLab id: {}'.format(self.fastq_folder))
                logging.debug('MyDemuxUnit: Infer project info: {}'.format(self.project))
            elif tmatches_I:
                self.project = tmatches_I.groups()[0]
                # some project does not have a iLab id, so leave it as '-1' while still keep this project
                logging.warning('MyDemuxUnit: Failed to infer iLab id: {}'.format(self.fastq_folder))
                logging.debug('MyDemuxUnit: Infer project info: {}'.format(self.project))
            else:
                # some project does not have a iLab id, so leave it as '-1' while still keep this project
                ###self.valid = False
                # get project name
                self.project = folder_name
                logging.warning('MyDemuxUnit: Failed to infer iLab id: {}'.format(self.fastq_folder))
                logging.debug('MyDemuxUnit: Infer project info: {}'.format(self.project))
            # save ilab as an integer
            self.ilab = int(self.ilab)

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
                logging.debug('MyDemuxUnit: Infer demux report html file: {}'.format(self.report_file))
            else:
                self.valid = False
                logging.warning('MyDemuxUnit: Failed to locate demux report html file: {}'.format(tsumfile))

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
            for root, directories, files in os.walk(self.fastq_folder):
                tfqs.extend([os.path.join(root, x) for x in files if search('fastq\.gz$',x) and not search('^Undetermined_',x)])
            # fastq files found?
            if len(tfqs) > 0:
                # read 1/2 fastq
                tr1s = [x for x in tfqs if search('R1_001\.fastq\.gz', x)]
                tr2s = [x for x in tfqs if search('R2_001\.fastq\.gz', x)]
                # find both read 1 and read 2 files
                if len(tr1s) == 0:
                    self.valid = False
                    logging.warning('MyDemuxUnit: read 1 file does not exist, failed to infer read length.')
                    return None
                # SR or PE
                if len(tr2s) == 0:
                    self.seqtype = 'SR'
                else:
                    self.seqtype = 'PE'
                # read length
                self.read_1_len = self.guess_read_length(tr1s[0])
                if len(tr2s) > 0:
                    self.read_2_len = self.guess_read_length(tr2s[0])
                else:
                    self.read_2_len = 0
                logging.debug('MyDemuxUnit: Infer sequencing type: {}'.format(self.seqtype))
                logging.debug('MyDemuxUnit: Infer read 1 length: {}'.format(self.read_1_len))
                if len(tr2s) > 0:
                    logging.debug('MyDemuxUnit: Infer read 2 length: {}'.format(self.read_2_len))

    def load_report(self):
        """extract demux summary from report file"""
        if self.report_file is not None:
            sumtables = pd.read_html(self.report_file, match='Sample')
            # double check on detected tables
            if len(sumtables) != 1:
                logging.warning("MyDemuxUnit: More than one table with 'Sample' column detected in demux summary file {}. Use the first table.".format(self.report_file))
            else:
                # use the first table
                mytable = sumtables[0]
                # merge reads from the same sample but different lanes
                self.summary = mytable.groupby(['Sample','Barcode sequence'])\
                    .agg(num_lanes=('PF Clusters','count'), PF_Clusters=('PF Clusters','sum'), Mbases=('Yield (Mbases)','sum'), Q30=('% >= Q30bases','mean'))\
                    .reset_index()\
                    .rename(columns={'PF_Clusters':'PF Clusters','Mbases':'Yield (Mbases)','Q30':'% >= Q30bases'})\
                #logging.debug('MyDemuxUnit: {}'.format(self.summary.head()))
            # double check if all required columns exist
            for col in self.columns:
                if col not in self.summary.columns:
                    self.valid = False
                    logging.warning("MyDemuxUnit: Column '{}' cannot be found in the demux summary file {}.".format(col, self.report_file))
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
        # a valid demux unit and no existing tables?
        if not self.valid:
            logging.warning('MyDemuxUnit: Not a valid demux project: {}'.format(self.fastq_folder))
            return None
        if self.detail is not None:
            logging.warning('MyDemuxUnit: existing detail table.')
            return None
        # get a copy of summary for output purpose
        self.detail = self.summary[self.columns].copy()
        # add other information
        self.detail.insert(0, 'iLab', [self.ilab] * self.detail.shape[0])
        self.detail.insert(1, 'project', [self.project] * self.detail.shape[0])
        self.detail.insert(2, 'seqtype', [self.seqtype] * self.detail.shape[0])
        self.detail.insert(3, 'read1len', [self.read_1_len] * self.detail.shape[0])
        self.detail.insert(4, 'read2len', [self.read_2_len] * self.detail.shape[0])
        # sort by iLab id
        self.detail.sort_values(by='iLab', inplace=True)

    def to_file(self, outfile):
        """combine information and write to file"""
        # detail table ready?
        if self.detail is None:
            return None
        # write to file
        if search('\.xlsx$', outfile):
            self.detail.to_excel(outfile, index=False)
        elif search('\.csv$', outfile):
            self.detail.to_csv(outfile, sep=',', index=False)
        elif search('\.tsv$|\.txt$', outfile):
            self.detail.to_csv(outfile, sep='\t', index=False)
        else:
            logging.warning('MyDemuxUnit: Unsupported output file extension: {}'.format(outfile.split('.')[-1]))
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

    def __init__(self, run_folder, date=None, flowcell_id=None, platform=None, seqtype=None, read_1_len=-1, read_2_len=-1, index_1_len=-1, index_2_len=-1, detail_table=None):
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
        # detailed data table
        self.detail = detail_table
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
            else:# otherwise, assume this is a novaseq run
                self.platform = 'novaseq'
            #    self.valid = False
            #    logging.warning('MyDemuxRun: Failed to infer platform for {}'.format(self.platform))
            logging.debug('MyDemuxRun: Infer platform: {}'.format(self.platform))

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
                logging.warning('MyDemuxRun: Failed to infer sequencing date for {}'.format(self.run_folder))
            else:
                self.date = '20{}-{}-{}'.format(tdate[:2],tdate[2:4],tdate[4:])
                logging.debug('MyDemuxRun: Infer sequencing date: {}'.format(self.date))

    def infer_flowcell_id(self):
        """guess flowcell id"""
        if self.fcid is None:
            self.fcid = os.path.basename(self.run_folder).split('_')[-1]
            if self.platform not in ['nextseq2000','miseq']:
                self.fcid = self.fcid[1:]# the first letter refers to flowcell A or B
            logging.debug('MyDemuxRun: Infer flowcell id: {}'.format(self.fcid))

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
            logging.warning('MyDemuxRun: Unable to read the run parameter file: {}'.format(run_para_file))
            return None
        # parse xml file
        tree = ET.parse(run_para_file)
        # get root
        root = tree.getroot()
        # extract read length info
        for x in root.iter(r1id):
            self.read_1_len = int(x.text)
            break
        for x in root.iter(r2id):
            self.read_2_len = int(x.text)
            break
        for x in root.iter(i1id):
            self.index_1_len = int(x.text)
            break
        for x in root.iter(i2id):
            self.index_2_len = int(x.text)
            break
        # do we have length for all reads?
        if self.read_1_len == -1:
            self.valid = False
            logging.warning('MyDemuxRun: Failed to get read 1 length.')
        elif self.read_2_len == -1:
            self.valid = False
            logging.warning('MyDemuxRun: Failed to get read 2 length.')
        elif self.index_1_len == -1:
            self.valid = False
            logging.warning('MyDemuxRun: Failed to get index 1 length.')
        elif self.index_2_len == -1:
            self.valid = False
            logging.warning('MyDemuxRun: Failed to get index 2 length.')
        else:
            logging.debug('MyDemuxRun: Infer read 1 length: {}'.format(self.read_1_len))
            logging.debug('MyDemuxRun: Infer read 2 length: {}'.format(self.read_2_len))
            logging.debug('MyDemuxRun: Infer index 1 length: {}'.format(self.index_1_len))
            logging.debug('MyDemuxRun: Infer index 2 length: {}'.format(self.index_2_len))

    def infer_seq_type(self):
        """infer sequencing type based on read length"""
        if self.seqtype is None:
            if self.read_2_len > 0:
                self.seqtype = 'PE'
            else:
                self.seqtype = 'SR'
            logging.debug('MyDemuxRun: Infer sequencing type: {}'.format(self.seqtype))

    def extract_demux_units(self):
        """extract demux information for each project"""
        # screen for fastq folders under run folder
        # assume fastq stored under folder named 'Unaligned', but not necessarily directly under it (e.g. 10X runs)
        candidate_folders = [x for x in os.listdir(self.run_folder) if search('Unaligned',x)]
        # scan for folders containing fastq files (demux unit)
        temp_folders = []
        for folder in candidate_folders:
            for root, directories, files in os.walk(os.path.join(self.run_folder, folder)):
                for tfile in files:
                    if search('\.fastq\.gz$', tfile) and not search('^Undetermined_', tfile):
                        # search for the sub-folder right underneath the 'UnalignedXX' or 'fastq_path' folder
                        # --- there are cases where fastq files were saved in a separate folder per sample
                        # --- for current runs, a fastq folder is named in the format 'PI-user-iLab(_20XX_XX_XX)'
                        # --- for older runs, the folder name may be arbitrary
                        # ---------------------------------------------------------------------------------------
                        # separate path into pieces
                        items = os.path.normpath(root).split(os.sep)
                        # search for 'fastq_path' or 'Unaligned'
                        k = len(items) - 1
                        while k >= 0:
                            if search('Unaligned|fastq_path', items[k]):
                                break
                            k -= 1
                        if k < 0:# failed to find Unaligned or fastq_path folder
                            logging.warning("MyDemuxRun: failed to locate fastq folder, unable to detect key words 'Unaligned' or 'fastq_path'")
                            # return the current folder as fastq_path, another downstream warning messages will give more details for debugging
                            temp_folders.append(root)
                        else:# find the key word
                            # use the folder right underneath it
                            temp_folders.append(os.sep.join(items[:k+2]))
                        break
        fastq_folders = list(set(temp_folders))
        #logging.debug('MyDemuxRun: '+'\n'.join(fastq_folders))
        # process each demux unit and store it as a MyDemuxUnit object
        for folder in fastq_folders:
            logging.debug('MyDemuxRun: {}'.format(folder))
            tunit = MyDemuxUnit(folder, self.fcid)
            tunit.infer_all()
            if tunit.valid:
                self.units.append(tunit)
        # at least one valid unit collected?
        if not self.units:
            self.valid = False
            logging.warning('MyDemuxRun: Not a valid demux run since no valid units (projects) are found: {}'.format(self.run_folder))

    def infer_all(self):
        """infer all available information"""
        self.infer_platform()
        self.infer_seq_date()
        self.infer_flowcell_id()
        self.infer_read_length()
        self.infer_seq_type()
        self.extract_demux_units()

    def dedup_records(self, mydf):
        """de-duplicate records within one demux run"""
        # sometimes, one project may be demuxed several times, each time with different demux settings.
        # this function performes de-duplication such that one record is kept per sample based on the following criteria:
        # 1) keep record with the highest number of reads
        # 2) for multiple sequencing types: choose SR if it applies
        # 3) for different read lengths, choose the short one

        # categorize the 'seqtype' column
        mydf['seqtype'] = pd.Categorical(mydf['seqtype'], ordered=True, categories=['SR','PE'])
        # sort by reads, sequencing type and read lengths
        mydf.sort_values(by=['PF Clusters','seqtype','read1len','read2len'], ascending=[False,True,True,True], inplace=True)
        # dedup by choosing the first record
        return mydf.drop_duplicates(subset=['iLab','Sample'], keep='first')

    def prepare_table(self):
        """combine information into a dataframe table"""
        # a valid demux run + no existing tables?
        if not self.valid or self.detail is not None:
            return None
        # merge information from each demux unit
        self.detail = pd.concat([x.prepare_table() for x in self.units])
        # add additional information
        # platform
        self.detail.insert(2, 'platform', [self.platform] * self.detail.shape[0])
        # flowcell id
        self.detail.insert(3, 'flowcell_id', [self.fcid] * self.detail.shape[0])
        # sequencing date
        self.detail.insert(4, 'date', [self.date] * self.detail.shape[0])
        self.detail['date'] = pd.to_datetime(self.detail['date'], format="%Y-%m-%d")
        # run sequencing type
        self.detail.insert(self.detail.shape[1], 'run_seqtype', [self.seqtype] * self.detail.shape[0])
        # run read length
        self.detail.insert(self.detail.shape[1], 'run_read1len', [self.read_1_len] * self.detail.shape[0])
        self.detail.insert(self.detail.shape[1], 'run_read2len', [self.read_2_len] * self.detail.shape[0])
        self.detail.insert(self.detail.shape[1], 'run_index1len', [self.index_1_len] * self.detail.shape[0])
        self.detail.insert(self.detail.shape[1], 'run_index2len', [self.index_2_len] * self.detail.shape[0])
        # de-dupliate record by sample
        self.detail = self.dedup_records(self.detail)
        # sort by sequencing date, instrument, iLab, project name, sample name
        self.detail.sort_values(by=['date','platform','iLab','project','Sample'], inplace=True)
        #logging.debug('MyDemuxRun: {}'.format(self.detail.shape))
        #logging.debug('MyDemuxRun: {}'.format(self.detail.head(2)))

    def to_file(self, outfile):
        """combine information and write to file"""
        # details table ready?
        if self.detail is None:
            return None
        # write to file
        if search('\.xlsx$', outfile):
            self.detail.to_excel(outfile, index=False)
        elif search('\.csv$', outfile):
            self.detail.to_csv(outfile, sep=',', index=False)
        elif search('\.tsv$|\.txt$', outfile):
            self.detail.to_csv(outfile, sep='\t', index=False)
        else:
            logging.warning('MyDemuxRun: Unsupported output file extension: {}'.format(outfile.split('.')[-1]))
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

    def __init__(self, server_folder, detail_table=None):
        """class constructor"""
        # server folder
        self.server_folder = server_folder
        # sequencing runs inside this folder
        self.runs = []
        # detailed data table
        self.detail = detail_table
        # a valid demux folder for output?
        self.valid = True

    def extract_demux_runs(self):
        """extract demux information for all available runs inside the folder"""
        logging.info('[Folder]\t{}'.format(self.server_folder))
        # screen for all available sequencing runs
        # e.g. 210901_A00814_0481_AHJ5T2DRXY
        run_folders = [x for x in os.listdir(self.server_folder) if search('\d+_[A-Z\d]+_\d+_[A-Za-z\d]+',x)]
        # process each demux run folder and store it as a MyDemuxRun object
        for folder in run_folders:
            logging.info('[Run]\t{}'.format(folder))
            trun = MyDemuxRun(os.path.join(self.server_folder,folder))
            trun.infer_all()
            if trun.valid:
                self.runs.append(trun)
            #logging.debug('MyDemuxFolder: {}'.format(trun))
        # at least one valid demux run?
        if not self.runs:
            self.valid = False
            logging.warning('MyDemuxFolder: Not a valid demux folder since no valid runs are found: {}'.format(self.server_folder))

    def prepare_table(self):
        """combine information into a dataframe table"""
        # a valid demux folder + no existing table?
        if not self.valid or self.detail is not None:
            return None
        # merge information from each demux run
        self.detail = pd.concat([x.prepare_table() for x in self.runs])
        # sort by sequencing date, instrument, iLab, project name, sample name
        self.detail.sort_values(by=['date','platform','iLab','project','Sample'], inplace=True)
        #logging.debug('MyDemuxFolder: {}'.format(self.detail.shape))
        #logging.debug('MyDemuxFolder: {}'.format(self.detail.head(2)))

    def to_file(self, outfile):
        """combine information and write to file"""
        # details table ready?
        if self.detail is None:
            return None
        # write to file
        if search('\.xlsx$', outfile):
            self.detail.to_excel(outfile, index=False)
        elif search('\.csv$', outfile):
            self.detail.to_csv(outfile, sep=',', index=False)
        elif search('\.tsv$|\.txt$', outfile):
            self.detail.to_csv(outfile, sep='\t', index=False)
        else:
            logging.warning('MyDemuxFolder: Unsupported output file extension: {}'.format(outfile.split('.')[-1]))
            return None
        logging.info('MyDemuxFolder: write MyDemuxFolder to file: {}'.format(outfile))

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

    def __init__(self, server_folder_list, detail_table=None, overview_table=None):
        """class constructor"""
        # a list of server folders to screen
        self.server_folder_list = server_folder_list
        # demux info for all available server folders
        self.folders = []
        # detailed data table
        self.detail = detail_table
        # overview data table
        self.overview = overview_table
        # columns to include
        self.columns = ['iLab','project','platform','flowcell_id','date','seqtype','read1len','read2len',\
            'Sample','Barcode sequence','num_lanes','PF Clusters','Yield (Mbases)','% >= Q30bases',\
            'run_seqtype','run_read1len','run_read2len','run_index1len','run_index2len']
        # month abbr.
        self.months = {1:'Jan',2:'Feb',3:'Mar',4:'Apr',5:'May',6:'Jun',7:'Jul',8:'Aug',9:'Sep',10:'Oct',11:'Nov',12:'Dec'}
        # a valid demux auto for output?
        self.valid = True

    def extract_demux_folders(self):
        """extract demux information for all available folders in the list"""
        for folder in self.server_folder_list:
            tfolder = MyDemuxFolder(folder)
            tfolder.extract_demux_runs()
            if tfolder.valid:
                self.folders.append(tfolder)
            #logging.debug('MyDemuxFolder: {}'.format(tfolder))
        # at least one valid demux folder?
        if not self.folders:
            self.valid = False
            logging.warning('MyDemuxFolder: Not a valid demux auto since no valid folders are found: \n{}'.format('\n'.join(self.server_folder_list)))

    def prepare_table(self):
        """combine information into a dataframe table"""
        # a valid demux auto + no existing tables?
        if not self.valid or self.detail is not None or self.overview is not None:
            return None
        # merge information from all demux folders
        self.detail = pd.concat([x.prepare_table() for x in self.folders])
        # sort by sequencing date, instrument, iLab, project name, sample name
        self.detail.sort_values(by=['date','platform','iLab','project','Sample'], inplace=True)
        # prepare an overview table at the project level
        self.overview = self.detail.groupby(by=['iLab','project','platform','flowcell_id','date'])['PF Clusters'].agg('sum').reset_index(name='PF Clusters')
        # sort the overview table by date and iLab
        self.overview.sort_values(by=['date','iLab'], inplace=True)
        #logging.debug('MyDemuxFolder: {}'.format(self.detail.shape))
        #logging.debug('MyDemuxFolder: {}'.format(self.detail.head(2)))

    def write_df_to_excel(self, mydf, outdir, outprefix):
        """write a DataFrame to an excel file, with sheets by months/overview"""
        # prepare an overview table at the project level
        overview = mydf.groupby(by=['iLab','project','platform','flowcell_id','date'])['PF Clusters'].agg('sum').reset_index(name='PF Clusters')
        # sort the overview table by date and iLab
        overview.sort_values(by=['date','iLab'], inplace=True)
        # output file
        outfile = os.path.join(outdir, '{}.{}.xlsx'.format(outprefix, mydf['year'].iloc[0]))
        # open an excel file handler and write to it
        with pd.ExcelWriter(outfile, datetime_format='mmm d, yyyy', engine='xlsxwriter') as writer:
            # get the xlsxwriter workbook
            workbook = writer.book
            # add some cell formats
            format_int = workbook.add_format({'num_format': '#,##'})
            # overview table
            overview.sort_values(by=['date','platform','iLab','project']).to_excel(writer, sheet_name='overview', index=False)
            # apply format to sheet 'overview'
            worksheet = writer.sheets['overview']
            worksheet.set_column('B:B', 20)
            worksheet.set_column('C:C', 12)
            worksheet.set_column('D:D', 15)
            worksheet.set_column('E:E', 13)
            worksheet.set_column('F:F', 15, format_int)
            worksheet.autofilter('A1:F{}'.format(overview.shape[0]+1))
            worksheet.freeze_panes(1, 0)
            # per-month information table
            mydf['month'] = mydf['date'].dt.month
            for m in range(12):
                select_month = (mydf['month'] == m+1)
                if select_month.sum() > 0:
                    sheet_name = self.months[m+1]
                    mydf[select_month][self.columns].sort_values(by=['date','platform','iLab','project','Sample']).to_excel(writer, sheet_name=sheet_name, index=False)
                    # apply format to current sheet
                    worksheet = writer.sheets[sheet_name]
                    worksheet.set_column('B:B', 15)
                    worksheet.set_column('C:C', 12)
                    worksheet.set_column('D:D', 15)
                    worksheet.set_column('E:E', 13)
                    worksheet.set_column('F:F', 12)
                    worksheet.set_column('G:G', 13)
                    worksheet.set_column('H:H', 13)
                    worksheet.set_column('I:I', 20)
                    worksheet.set_column('J:J', 22)
                    worksheet.set_column('K:K', 15)
                    worksheet.set_column('L:L', 15, format_int)
                    worksheet.set_column('M:M', 18, format_int)
                    worksheet.set_column('N:N', 18)
                    worksheet.set_column('O:O', 16)
                    worksheet.set_column('P:P', 16)
                    worksheet.set_column('Q:Q', 16)
                    worksheet.set_column('R:R', 16)
                    worksheet.set_column('S:S', 16)
                    worksheet.autofilter('A1:S{}'.format(select_month.sum()+1))
                    worksheet.freeze_panes(1, 0)
            # Close the Pandas Excel writer and output the Excel file.
            writer.save()

    def to_excel(self, outdir, outprefix='GRCF.demux.summary'):
        """combine information and write to an excel file"""
        # unlike writting to a plain text file
        # 1) one excel per year, one spread sheet per month
        # 2) include an overview spread sheet storing the total reads per project

        # overview/details tables ready?
        if self.detail is None or self.overview is None:
            return None
        # separate data by year and write to file
        self.detail['year'] = self.detail['date'].dt.year
        # sort data by year
        #self.detail.sort_values(by=['year'], inplace=True)
        self.detail.groupby(by='year').apply(self.write_df_to_excel, outdir=outdir, outprefix=outprefix)

    def to_file(self, outdir, outprefix='GRCF.demux.summary', outext='txt'):
        """combine information and write to file"""
        # overview/details tables ready?
        if self.detail is None or self.overview is None:
            return None
        # write to file
        if outext == 'xlsx':
            self.to_excel(outdir, outprefix)
            logging.info('write MyDemuxAuto to excel files: {}'.format(os.path.join(outdir, '{}.XXXX.xlsx'.format(outprefix))))
        else:
            detail_file = os.path.join(outdir, '{}.details.{}'.format(outprefix, outext))
            overview_file = os.path.join(outdir, '{}.overview.{}'.format(outprefix, outext))
            if outext == 'csv':
                self.detail.to_csv(detail_file, sep=',', index=False)
                self.overview.to_csv(overview_file, sep=',', index=False)
            elif outext == 'tsv' or outext == 'txt':
                self.detail.to_csv(detail_file, sep='\t', index=False)
                self.overview.to_csv(overview_file, sep='\t', index=False)
            else:
                logging.warning('MyDemuxFolder: Unsupported output file extension: {}'.format(outfile.split('.')[-1]))
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
    #du = MyDemuxUnit('/data/seq/NovaSeq6000/200625_A00814_0207_BHMNKGDRXX/Unaligned_1/Landau-RC-8951_2020_06_25','HMNKGDRXX')
    #du.infer_ilab()
    #du.infer_report_file()
    #du.infer_seq_type()
    #du.load_report()
    du.infer_all()
    du.prepare_table()
    du.to_file('/tmp/test.txt')
    print(du)
    #print(du.summary.columns)

def test_MyDemuxRun():
    #dr = MyDemuxRun('/gc7-data/NovaSeq6000/211015_A00814_0503_AHL5G2DSX2')
    #dr = MyDemuxRun('/scratch/seq_data/NovaSeq6000/211014_A00814_0502_BHL5H3DSX2')
    #dr = MyDemuxRun('/data/seq/NovaSeq6000/200625_A00814_0207_BHMNKGDRXX')
    #dr = MyDemuxRun('/gc7-data/NovaSeq6000/210730_A00814_0465_AHHKGVDSX2')
    dr = MyDemuxRun('/gc-archive2/gc5-backup/GRCF_data_archive/NovaSeq6000/200309_A00814_0164_AHNMV7DMXX')
    #dr = MyDemuxRun('/gc-archive2/gc5-backup/GRCF_data_archive/NovaSeq6000/200724_A00814_0224_BHMHTLDRXX')
    #dr.infer_platform()
    #dr.infer_seq_date()
    #dr.infer_flowcell_id()
    #dr.infer_read_length()
    #dr.infer_seq_type()
    #dr.extract_demux_units()
    #dr.prepare_table()
    dr.infer_all()
    dr.prepare_table()
    dr.to_file('/tmp/test.run.txt')
    print(dr)
    #print(len(dr.units))
    #print(dr.units[0])

def test_MyDemuxFolder():
    #df = MyDemuxFolder('/scratch/seq_data/NovaSeq6000')
    df = MyDemuxFolder('/gc7-data/NovaSeq6000')
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
    da.to_file('/tmp','test.auto','xlsx')
    #da.to_file('/tmp','test.auto','txt')
    #da.to_excel('/tmp','test.auto')

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

def run_MyDemuxAuto():
    """scan GRCF server demux folders"""
    # run on the gc6 server
    server_folders = ['/data/seq/NovaSeq6000','/scratch/seq_data/NovaSeq6000','/data/seq/gc7_demux_runs',\
        '/gc7-data/NovaSeq6000','/gc7-data/NextSeq500',\
        '/gc5/NextSeq500','/gc5/NovaSeq6000',\
        '/gc4/NextSeq2000/NextSeq2000','/gc4/NextSeq500','/gc4/HiSeq2500_new/flowcellA','/gc4/HiSeq2500_new/flowcellB',\
        '/gc-archive2/gc4-backup/HiSeq4000/flowcellA','/gc-archive2/gc4-backup/GRCF_data_archive/HiSeq4000',\
        '/gc-archive2/gc4-backup/GRCF_data_archive/NextSeq2000','/gc-archive2/gc4-backup/GRCF_data_archive/NextSeq500',\
        '/gc-archive2/gc5-backup/GRCF_data_archive/NovaSeq6000','/genome2/GRCF_data_archive/HiSeq2500',\
        '/genome2/GRCF_data_archive/HiSeq4000','/genome2/GRCF_data_archive/NextSeq500']
    da = MyDemuxAuto(server_folders)
    da.extract_demux_folders()
    da.prepare_table()
    da.to_file('/data/seq/tmp','GRCF.demux.summary','xlsx')
    da.to_file('/data/seq/tmp','GRCF.demux.summary','txt')

def main():
    # set up logging
    root_logger = setup_logging('/data/seq/tmp/GRCF.demux.summary.auto.log', level=logging.INFO)
    #root_logger = setup_logging('/tmp/test.auto.log', level=logging.DEBUG)
    #root_logger = setup_logging('/tmp/test.auto.log', level=logging.INFO)

    #test_MyDemuxUnit()
    #test_MyDemuxRun()
    #test_MyDemuxFolder()
    #test_MyDemuxAuto()
    run_MyDemuxAuto()

# main
if __name__ == '__main__':
    main()

