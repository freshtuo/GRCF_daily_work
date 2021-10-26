#!/usr/bin/env python
# collect_demux_info.py
# collect demux info from sequencing runs, for records
#

import sys
import gzip
import pandas as pd

from os.path import exists
from os import getcwd
from os import listdir
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

    def __init__(self, fastq_folder, flowcell_id, report_file=None, ilab=None, project=None, seqtype=None, read_1_len=-1, read_2_len=-1):
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

    def infer_ilab(self):
        """guest iLab id and project name"""
        # eg: Loda-MJ-10557_2021_06_15 or Diaz-Meco-MADM-10830_2021_07_19
        if self.ilab is None:
            folder_name = self.fastq_folder.split('/')[-1]
            tmatches = findall('(\d+)', folder_name)
            if tmatches:
                self.ilab = int(tmatches[0])
                self.project = folder_name[:folder_name.find(tmatches[0])-1]
                print('Infer iLab id: {}'.format(self.ilab))
                print('Infer PI info: {}'.format(self.project))

    def infer_report_file(self):
        """guess demux report html file if not provided"""
        # report file already exists?
        if self.report_file is None:
            # the 'Reports' folder should be in the same level as the fastq folder and 
            # they share the same parent folder, not necessarily named 'Unaligned' (e.g. 10X runs)
            tfolders = self.fastq_folder.split('/')
            # guess report html file
            tsumfile = '/'.join(tfolders[:-1]+['Reports','html',self.fcid,tfolders[-1].split('_')[0],'all','all','laneBarcode.html'])
            # check file existence
            if exists(tsumfile):
                self.report_file = tsumfile
                print('Infer demux report html file: {}'.format(self.report_file))

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
            tfqs = [x for x in listdir(self.fastq_folder) if search('fastq\.gz',x)]
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
                self.read_1_len = self.guess_read_length('{}/{}'.format(self.fastq_folder, tr1s[0]))
                if len(tr2s) > 0:
                    self.read_2_len = self.guess_read_length('{}/{}'.format(self.fastq_folder, tr2s[0]))
                print('Infer sequencing type: {}'.format(self.seqtype))
                print('Infer read 1 length: {}'.format(self.read_1_len))
                if len(tr2s) > 0:
                    print('Infer read 2 length: {}'.format(self.read_2_len))

    def load_report(self):
        """extract demux summary from report file"""
        if self.report_file is not None:
            sumtables = pd.read_html(self.report_file, match='Sample')
            # double check on detected tables
            if len(sumtables) != 1:
                print("Warning: more than one table with 'Sample' column detected in demux summary file {}".format(self.report_file))
            else:
                self.summary = sumtables[0]

    def to_file(self, outfile):
        """combine information and write to file"""
        # get a copy of summary for output purpose
        mytable = self.summary[self.columns].copy()
        # add other information
        mytable.insert(0, 'iLab', [self.ilab] * mytable.shape[0])
        mytable.insert(1, 'project', [self.project] * mytable.shape[0])
        mytable.insert(2, 'seqtype', [self.seqtype] * mytable.shape[0])
        mytable.insert(3, 'read1len', [self.read_1_len] * mytable.shape[0])
        mytable.insert(4, 'read2len', [self.read_2_len] * mytable.shape[0])
        # write to file
        if search('\.xlsx$', outfile):
            mytable.to_excel(outfile, index=False)
        elif search('\.csv$', outfile):
            mytable.to_csv(outfile, sep=',', index=False)
        elif search('\.tsv$|\.txt$', outfile):
            mytable.to_csv(outfile, sep='\t', index=False)
        else:
            print('Unsupported output file extension: {}'.format(outfile.split('.')[-1]))

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
        return '\n'.join(to_print)

# MyDemuxRun: a sequencing run
#class MyDemuxRun:
    

# functions
def main():
    du = MyDemuxUnit('/gc7-data/NovaSeq6000/211015_A00814_0503_AHL5G2DSX2/Unaligned_1/Guzman-NMT-11319_2021_10_15','HL5G2DSX2')
    du.infer_ilab()
    du.infer_report_file()
    du.infer_seq_type()
    du.load_report()
    du.to_file('/tmp/test.txt')
    print(du)
    #print(du.summary.columns)

# main
if __name__ == '__main__':
    main()

