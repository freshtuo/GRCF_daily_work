#!/usr/bin/env python
# review_samplesheet_index.py
# examine the index sequence conflicts in a samplesheet file
# 

import sys
import pandas as pd
import os.path
import logging

from os.path import exists
from re import search
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter

# class
class MySampleSheet:
    """Class for collecting data from a samplesheet file"""

    def __init__(self, args):
        """class constructor"""
        # samplesheet filename
        self.filename = args.samplesheetfile
        # date
        self.date = None
        # samplesheet header lines
        self.header = []
        # samplesheet table
        self.table = None
        # original columns in samplesheet
        self.orig_columns = []
        # 10X ref index folders (/data/seq/QC_Share/GRCF/Sequencing/SampleSheets/10X_indexes_all)
        self.tenx_dir = args.tenx_dir
        # 10X ref index file
        self.tenx_idx_file = args.tenx_idx_file
        # 10X index
        self.tenx_idx = None
        # fixed samplesheet table
        self.fixed_table = None
        # columns to be considered as string
        self.columns_as_str = ['Sample_ID','Sample_Name','I7_Index_ID','index','I5_Index_ID','index2','Sample_Project']
        # names for index columns
        self.i7_index_id = 'I7_Index_ID'
        self.i7_index_seq = 'index'
        self.i5_index_id = 'I5_Index_ID'
        self.i5_index_seq = 'index2'
        # column to search for 10X index ids
        self.si_index_id = 'index'
        # output details?
        self.detail = args.detail
        # updated samplesheet file
        self.updatedsheet = args.updatedsheetfile

    def get_info(self):
        """extract lanes and samples information"""
        logging.info('Extracting lanes and samples information.')
        # file exists?
        if not exists(self.filename):
            logging.error('failed.\nUnable to read the samplesheet file: {}'.format(self.filename))
            sys.exit(1)
        # extract header lines
        k = 0
        with open(self.filename, 'r') as fin:
            # read entire sheet
            data = fin.readlines()
            # locate the start line of DATA section
            while k < len(data):
                if search('^\[Data\]', data[k]):
                    break
                k += 1
            # failed to locate DATA section
            if k == len(data):
                logging.error('failed.\nFailed to locate [DATA] section in samplesheet file: {}'.format(self.filename))
                sys.exit(2)
            # save header lines
            self.header.extend(data[:k+1])
        # read sample table as pandas DataFrame
        self.table = pd.read_table(self.filename, sep=',', header=0, skiprows=k+1, low_memory=False)
        # replace NAs with ''
        self.table.fillna('', inplace=True)
        # record original columns (for output purpose)
        self.orig_columns.extend(self.table.columns.tolist())
        # add 'Lane' column if it does not exists, and fill in 1
        if 'Lane' not in self.table.columns:
            self.table['Lane'] = 1
        # convert columns to str in case any of them include numbers only
        for column in self.columns_as_str:
            self.table[column] = self.table[column].astype(str)

    def diff_two_indexes(self, idx1, idx2):
        """calculate Hamming distance between two indexes"""
        return len([i for i in range(len(idx1)) if idx1[i] != idx2[i]])

    def diff_list_indexes(self, idxList):
        """calculate pairwise Hamming distance on a list of indexes"""
        distances = []
        for i in range(len(idxList)):
            for j in range(i+1, len(idxList)):
                distances.append((i+1,j+1,idxList[i],idxList[j],self.diff_two_indexes(idxList[i],idxList[j])))
        return distances

    def load_10X_single_index(self, filename):
        """Load a 10X Single Index kit"""
        # read in csv file
        idx = pd.read_table(os.path.join(self.tenx_dir, filename), sep=',', header=None, names=[self.i7_index_id,'indexA','indexB','indexC','indexD'], low_memory=False)
        # reformat four indexes into multiple entries
        temp_list = []
        def myformat(x, temp_list):
            temp_list.append(pd.DataFrame({self.i7_index_id:[x[self.i7_index_id]]*4, 'index':[x['indexA'],x['indexB'],x['indexC'],x['indexD']], 'index2':['']*4}))
        idx.apply(func=myformat, axis=1, temp_list=temp_list)
        # combine indexes
        idx2 = pd.concat(temp_list)
        if self.detail:
            logging.info('{} barcodes from {} index ids loaded from {}.'.format(idx.shape[0], idx2.shape[0], filename))
        # return
        return idx2

    def load_10X_dual_index(self, filename):
        """Load a 10X Dual Index kit"""
        # read in csv file
        idx = pd.read_table(os.path.join(self.tenx_dir, filename), sep=',', header=0, comment='#', low_memory=False)
        # ignore i5 index for workflow A
        idx.drop(columns=['index2_workflow_a(i5)'], inplace=True)
        # rename columns for consistency
        idx.rename(columns={'index_name':self.i7_index_id,'index(i7)':'index','index2_workflow_b(i5)':'index2'}, inplace=True)
        if self.detail:
            logging.info('{} barcodes from {} index ids loaded from {}.'.format(idx.shape[0], idx.shape[0], filename))
        # return
        return idx

    def load_10X_index(self):
        """load 10X index (id + sequences)"""
        if self.tenx_idx_file is None and self.tenx_dir is None:
            logging.error('Unable to find 10X reference index sequences.')
            sys.exit(2)
        elif self.tenx_dir is not None:
            logging.info('Loading 10X reference index sequences from folder: {}.'.format(self.tenx_dir))
            idx_list = []
            # Dual Index kit NN
            idx_list.append(self.load_10X_dual_index('Dual_Index_Kit_NN_Set_A.csv'))
            # Dual Index kit NT
            idx_list.append(self.load_10X_dual_index('Dual_Index_Kit_NT_Set_A.csv'))
            # Dual Index kit TN
            idx_list.append(self.load_10X_dual_index('Dual_Index_Kit_TN_Set_A.csv'))
            # Dual Index kit TS
            idx_list.append(self.load_10X_dual_index('Dual_Index_Kit_TS_Set_A.csv'))
            # Dual Index kit TT
            idx_list.append(self.load_10X_dual_index('Dual_Index_Kit_TT_Set_A.csv'))
            # Single Index kit N
            idx_list.append(self.load_10X_single_index('Single_Index_Kit_N_Set_A.csv'))
            # Single Index kit T
            idx_list.append(self.load_10X_single_index('Single_Index_Kit_T_Set_A.csv'))
            # combine all indexes
            self.tenx_idx = pd.concat(idx_list)
            # if both 10X folder and index file are provided,
            # load index sequences from folder and overwrite the index file
            if self.tenx_idx_file is not None:
                logging.info('Write 10X reference index sequences to file: {}.'.format(self.tenx_idx_file))
                self.tenx_idx.to_csv(self.tenx_idx_file, sep='\t', index=False)
        else:
            logging.info('Loading 10X reference index sequences from file: {}.'.format(self.tenx_idx_file))
            self.tenx_idx = pd.read_table(self.tenx_idx_file, header=0, sep='\t', low_memory=False)
            # fill na entries
            self.tenx_idx.fillna({'index':'','index2':''}, inplace=True)

    def update_index(self):
        """update the missing 10X index in the samplesheet"""
        logging.info('Filling up missing 10X index sequences.')
        # divide sample sheet into two parts: with and without the need to fix index sequences
        data_to_fix = self.table[self.table[self.si_index_id].str.contains('^SI')]
        data_to_keep = self.table[~self.table[self.si_index_id].str.contains('^SI')]
        # copy index id to the right column
        data_to_fix[self.i7_index_id] = data_to_fix[self.si_index_id]
        # fill in index sequences
        data_to_fix = data_to_fix.drop(columns=['index','index2']).merge(self.tenx_idx, how='left', on=self.i7_index_id)
        # merge two parts
        self.fixed_table = pd.concat([data_to_keep, data_to_fix[data_to_keep.columns]])

    def screen_duplicate_index(self, td):
        """search for duplicate index"""
        # combine I7 and I5
        td['combined_index'] = td.apply(lambda x: '{}-{}'.format(x[self.i7_index_seq], x[self.i5_index_seq]), axis=1)
        # any identical indexes?
        duplicate_entries = td['combined_index'].duplicated(keep=False)
        if duplicate_entries.sum() > 0:
            logging.error('Duplicate indexes detected in lane {}:'.format(td['Lane'].iloc[0]))
            logging.error(td[duplicate_entries])

    def trim_index(self, idx_list):
        """trim a list of indexes to the same length"""
        # discard empty indexes
        idx_list = [x for x in idx_list if x != '']
        # minimum index length
        if idx_list:
            min_len = min([len(x) for x in idx_list])
            return [x[:min_len] for x in idx_list]
        else:
            return []

    def screen_close_index(self, td):
        """search for indexes with similar sequences (Hamming distance < 3)"""
        # i7 index
        # pairwise distance
        i7_dist = self.diff_list_indexes(self.trim_index(td[self.i7_index_seq].drop_duplicates(keep='first').tolist()))
        # minimum distance to return
        i7_dist_min = 99
        if i7_dist:
            i7_dist_min = min([x[-1] for x in i7_dist])
        # print out identical indexes in any
        i7_dist_same = [x for x in i7_dist if x[-1] == 0]
        if i7_dist_same:
            logging.error('identical i7 indexes found in lane {} after trimming bases to the same length.'.format(td['Lane'].iloc[0]))
            logging.error('identical i7 indexes:\n'+'\n'.join(['idxA: {}; idxB: {}; distance: {}'.format(x[2],x[3],x[4]) for x in i7_dist_same]))
        # print out close indexes if any
        i7_dist_close = [x for x in i7_dist if x[-1] < 3]
        if i7_dist_close:
            logging.warning('minimum distance for i7 index is smaller than 3 for lane {}, no mismatches allowed.'.format(td['Lane'].iloc[0]))
            if self.detail:
                logging.warning('close indexes:\n'+'\n'.join(['idxA: {}; idxB: {}; distance: {}'.format(x[2],x[3],x[4]) for x in i7_dist_close]))
        # i5 index
        # pairwise distance
        i5_dist = self.diff_list_indexes(self.trim_index(td[self.i5_index_seq].drop_duplicates(keep='first').tolist()))
        # minimum distance to return
        i5_dist_min = 99
        if i5_dist:
            i5_dist_min = min([x[-1] for x in i5_dist])
        # print out identical indexes in any
        i5_dist_same = [x for x in i5_dist if x[-1] == 0]
        if i5_dist_same:
            logging.error('identical i5 indexes found in lane {} after trimming bases to the same length.'.format(td['Lane'].iloc[0]))
            logging.error('identical i5 indexes:\n'+'\n'.join(['idxA: {}; idxB: {}; distance: {}'.format(x[2],x[3],x[4]) for x in i5_dist_same]))
        # print out close indexes if any
        i5_dist_close = [x for x in i5_dist if x[-1] < 3]
        if i5_dist_close:
            logging.warning('minimum distance for i5 index is smaller than 3 for lane {}, no mismatches allowed.'.format(td['Lane'].iloc[0]))
            if self.detail:
                logging.warning('close indexes:\n'+'\n'.join(['idxA: {}; idxB: {}; distance: {}'.format(x[2],x[3],x[4]) for x in i5_dist_close]))
        # return minimum distance
        return pd.Series([i7_dist_min, i5_dist_min], index=['i7','i5'])

    def examine_index(self):
        """examine whether or not there are conflicting indexes in the samplesheet"""
        logging.info('Examing conflicting indexes.')
        # 1. any identical indexes after combining I7 and I5 indexes
        self.fixed_table.groupby(by='Lane').apply(self.screen_duplicate_index)
        # 2. Hamming distance < 3 for either I7 or I5 indexes
        min_dist_per_lane = self.fixed_table.groupby(by='Lane').apply(self.screen_close_index)
        if self.detail:
            logging.info('Minimum index difference for i7 and i5 in each lane:\n{}'.format(min_dist_per_lane))

    def write_to_file(self):
        """write the fixed samplesheet to file"""
        if self.fixed_table is not None and self.updatedsheet is not None:
            with open(self.updatedsheet, 'w') as fsp:
                fsp.write(''.join(self.header))
                self.fixed_table.to_csv(fsp, sep=',', columns=self.orig_columns, index=False)

# functions
def get_arguments():
    """fetch command line arguments."""
    parser = ArgumentParser(description="""Given a samplesheet, check index sequence conflicts.""",
                            prog='review_samplesheet_index.py')
    parser.add_argument("-v", "--version", action="version", version='%(prog)s v0.1')
    parser.add_argument("-i", "--samplesheet", nargs="?", required=True, help="samplesheet file", metavar="samplesheet_file", dest="samplesheetfile")
    parser.add_argument("-l", "--log", nargs="?", default="review_samplesheet_index.log", help="log file", metavar="log_file", dest="logfile")
    parser.add_argument("-d", "--detail", action='store_true', default=False, help="output details", dest="detail")
    parser.add_argument("-o", "--updatedsheet", nargs="?", default=None, help="updated samplesheet file for manual check", metavar="updated_samplesheet_file", dest="updatedsheetfile")
    parser.add_argument("-t", "--tenxdir", nargs="?", default=None, help="load 10X index sequences in the given folder", metavar="tenx_dir", dest="tenx_dir")
    parser.add_argument("-r", "--tenxidx", nargs="?", default=None, help="load 10X index sequences from the given file", metavar="tenx_idx_file", dest="tenx_idx_file")
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
    # create a MySampleSheet object
    d = MySampleSheet(args)
    # load 10X index white list
    d.load_10X_index()
    # extract information from samplesheet
    d.get_info()
    # replace missing 10X indexes
    d.update_index()
    # write updated samplesheet to file
    d.write_to_file()
    # examing indexes
    d.examine_index()

# main
if __name__ == '__main__':
    main()

