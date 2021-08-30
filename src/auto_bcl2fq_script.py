#!/usr/bin/env python
# auto_bcl2fq_script.py
# automatically generate script for running demux given a samplesheet file
# by splitting samples based on index and sequencing type
# 

import sys
import pandas as pd

import xml.etree.ElementTree as ET

from os.path import exists
from os import getcwd
from re import search
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter

# class
class MyDemux:
    """Class for preparing demux"""

    def __init__(self, args):
        """class constructor"""
        # samplesheet header lines
        self.header = []
        # samplesheet table
        self.table = None
        # samplesheet mismatch allowance
        self.mismatch = None
        # read length
        self.read_1_len = -1
        self.read_2_len = -1
        self.index_1_len = -1
        self.index_2_len = -1
        # run folder
        self.run_folder = args.runfolder
        # input samplesheet file
        self.samplesheet_file = args.infile
        # skiplanes
        self.skip_lanes = args.skip
        # bcl2fastq program
        self.bcl2fastq = args.bcl2fastq
        # output script folder
        if args.scriptdir is None:
            self.script_folder = self.run_folder
        else:
            self.script_folder = args.scriptdir
        # output samplesheet folder
        if args.samplesheetdir is None:
            self.samplesheet_folder = '{}/Data/Intensities/BaseCalls'.format(self.run_folder)
        else:
            self.samplesheet_folder = args.samplesheetdir
        # overwrite existing output script file?
        self.overwrite = args.force
        # no lane splitting?
        self.no_lane_split = args.nosplit
        # platform
        self.platform = args.platform
        ## illegal characters (regular expression)
        #self.ilch = """[\?\(\)\[\]/\\=+<>:;"'\,\*\^|&\. ]"""
        # unusual characters (regular expression)
        self.nmch = '[^A-Z,a-z,0-9,\_,\-]'
        # original columns in the samplesheet sample table
        self.orig_columns = []
        # columns to check illegal characters
        self.columns_to_check = ['Sample_ID','Sample_Name','index','index2','Sample_Project','Description']

    def get_read_length(self):
        """extract sequenced read length from runparameter.xml"""
        print('Extracting read length...', end='')
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
        if not exists(run_para_file):
            print('failed.\nUnable to read the run parameter file: {}'.format(run_para_file))
            sys.exit(7)
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
            print('failed.\nFailed to get read 1 length.')
            sys.exit(8)
        if self.read_2_len == -1:
            print('failed.\nFailed to get read 2 length.')
            sys.exit(9)
        if self.index_1_len == -1:
            print('failed.\nFailed to get index 1 length.')
            sys.exit(10)
        if self.index_2_len == -1:
            print('failed.\nFailed to get index 2 length.')
            sys.exit(11)
        print('ok.')

    def show_read_length(self):
        print('read 1: {}'.format(self.read_1_len))
        print('read 2: {}'.format(self.read_2_len))
        print('index 1: {}'.format(self.index_1_len))
        print('index 2: {}'.format(self.index_2_len))

    def load_samplesheet(self):
        """load data from samplesheet"""
        print('Loading samplesheet...', end='')
        # samplesheet file exists?
        if not exists(self.samplesheet_file):
            print('failed.\nUnable to read the input samplesheet file: {}'.format(self.samplesheet_file))
            sys.exit(1)
        # extract header lines
        k = 0
        with open(self.samplesheet_file, 'r') as fin:
            # read entire sheet
            data = fin.readlines()
            # locate the start line of DATA section
            while k < len(data):
                if search('^\[Data\]', data[k]):
                    break
                k += 1
            # failed to locate DATA section
            if k == len(data):
                print('failed.\nFailed to locate [DATA] section in samplesheet file: {}'.format(self.samplesheet_file))
                sys.exit(2)
            # save header lines
            self.header.extend(data[:k+1])
        # read sample table as pandas DataFrame
        self.table = pd.read_table(self.samplesheet_file, sep=',', header=0, skiprows=k+1, low_memory=False)
        # replace NAs with ''
        self.table.fillna('', inplace=True)
        # record original columns (for output purpose)
        self.orig_columns.extend(self.table.columns.tolist())
        # add 'Lane' column if it does not exists, and fill in 1
        if 'Lane' not in self.table.columns:
            self.table['Lane'] = 1
        # skip lanes as desired
        self.table = self.table[~self.table['Lane'].isin(self.skip_lanes)]
        print('ok.')

    def column_duplicated(self, group, column):
        """check for duplicated items in a given column per group, in the samplesheet sample table"""
        # count duplicated items in the given column per given group
        # calculate overall counts
        ##dup_counts = self.table.groupby(by=group).apply(lambda x: x.duplicated(column).sum()).reset_index(name='counts')
        # label each entry according to whether or not it's a duplicate
        dup_marks = self.table.groupby(by=group).apply(lambda x: x.duplicated(column)).rename_axis(group+['index']).reset_index(name='duplicate')#.set_index('index')
        # subset duplicated entries
        dups = self.table[group+column][self.table.index.isin(dup_marks['index'][dup_marks['duplicate']].tolist())]
        # any duplicates?
        if not dups.empty:
            print('failed.\nDuplicated items found in columns {} for {}:'.format(' + '.join(column),' + '.join(group)))
            print(dups)
            sys.exit(4)

    def index_missing(self, df):
        """check if there's any index missing in a lane (for 'groupby.apply')"""
        # combine two indexes
        combined_index = df['index']+df['index2']
        # more than one sample (entry) with no indexes
        no_index = combined_index == ''
        if no_index.sum() > 1:
            return True
        # mixing samples of no indexes and samples of indexes
        if no_index.sum() > 0 and combined_index.size - no_index.sum() > 0:
            return True
        return False

    def examine_samplesheet(self):
        """examine samplesheet and make sure no illegal characters or duplicate items in it"""
        print('Examining samplesheet...', end='')
        # any illegal characters in selected columns
        for column in self.columns_to_check:
            tpats = self.table[column].str.contains(self.nmch)
            if tpats.sum() > 0:# unusual characters found!
                print('failed.\nUnusual characters found in column {}'.format(column))
                print(self.table[column][tpats])
                sys.exit(3)
        # check duplicate sample ids?
        self.column_duplicated(['Lane'], ['Sample_ID'])
        # check if indexes are unique?
        self.column_duplicated(['Lane'], ['index','index2'])
        # check if a project ('Sample_Project') has multiple sequencing types ('Description')
        tdes = self.table.drop_duplicates(['Sample_Project','Description']).groupby('Sample_Project')['Description'].count().reset_index(name='counts')
        tdes = tdes[tdes['counts'] > 1]
        if not tdes.empty:
            print('failed.\nMultiple sequencing types found in the following projects:')
            print(tdes['Sample_Project'])
            sys.exit(5)
        # check if indexes are missing?
        # no index is only allowed if there's only one sample in a lane
        missing_index = self.table.groupby('Lane').apply(self.index_missing).reset_index(name='Miss')
        missing_index = missing_index[missing_index['Miss']]
        if not missing_index.empty:
            print('failed.\nMissing index entries found in the following lanes:')
            print(missing_index)
            sys.exit(6)
        print('ok.')

    def diff(self,idx1, idx2):
        """calculate Hamming distance between two indexes"""
        return len([i for i in range(len(idx1)) if idx1[i] != idx2[i]])

    def minDiff(self,idxList):
        """calculate the minimum distance within a set of indexes"""
        # only one index in the given set? then return the length of that index
        if len(idxList) == 1:
            return len(idxList[0])
        # more than one index, do pairwise comparison and choose the minimum
        distances = []
        for i in range(len(idxList)):
            for j in range(i+1, len(idxList)):
                distances.append(self.diff(idxList[i], idxList[j]))
        return min(distances)

    def calculate_mismatch(self):
        """calculate the allowed mismatch for each Lane + index + project"""
        # rules:
        # 1) samples of same index type from the same project should be assigned the same mismatch allowance; 
        #    otherwise fastq data would be split into different folders (demux runs).
        # 2) do NOT allow mismatches for lanes with multiple index types
        print('Calculate allowed mismatches...', end='')
        # get index length
        self.table['index_len'] = self.table['index'].str.len()
        self.table['index2_len'] = self.table['index2'].str.len()
        # combine the two indexes
        self.table['combined_index'] = self.table.agg(lambda x:'{}{}'.format(x['index'],x['index2']), axis=1)
        # group entries by lane + index length
        # then calculate the minimum index distance per group
        mindists = self.table.groupby(by=['Lane','index_len','index2_len']).apply(lambda x: self.minDiff(x['combined_index'].tolist())).reset_index(name='mindist')
        # assign allowed mismatches
        # 1) initialize allowed mismatch based on minimum distance
        mindists['mismatch_1'] = mindists.apply(lambda x:x['mindist'] > 3, axis=1)
        # 2) multiple index type in a lane? then no mismatches
        mindists = mindists.merge(mindists.groupby(by=['Lane'])['Lane'].count().reset_index(name='num_index_type'), how='left', on='Lane')
        mindists['mismatch_2'] = mindists.apply(lambda x:x['num_index_type'] == 1, axis=1)
        # 3) use consistent allowed mismatches for the same project + index type
        # conclusion should be drawn based on mismatch_1 & mismatch_2
        mindists['mismatch_1&2'] = mindists['mismatch_1'] & mindists['mismatch_2']
        # add Sample_Project info (one column)
        mindists = self.table[['Lane','index_len','index2_len','Sample_Project']].drop_duplicates(['Lane','index_len','index2_len','Sample_Project']).merge(mindists, how='left', on=['Lane','index_len','index2_len'])
        mindists = mindists.merge(mindists.groupby(by=['index_len','index2_len','Sample_Project'])['mismatch_1&2'].all().reset_index(name='mismatch_3'), how='left', on=['index_len','index2_len','Sample_Project'])
        # final assignment by taking and on all three filters
        mindists['allow_mismatches'] = mindists['mismatch_1'] & mindists['mismatch_2'] & mindists['mismatch_3']
        # add to the sample table
        self.table = self.table.merge(mindists[['Lane','index_len','index2_len','Sample_Project','allow_mismatches']], how='left', on=['Lane','index_len','index2_len','Sample_Project'])
        # add the originally requested read type/length
        mindists = mindists.merge(self.table[['Lane','index_len','index2_len','Sample_Project','allow_mismatches','Description']].drop_duplicates(), how='left', on=['Lane','index_len','index2_len','Sample_Project','allow_mismatches'])
        # save mismatch filters
        self.mismatch = mindists[['Lane','index_len','index2_len','Sample_Project','mindist','num_index_type','mismatch_1','mismatch_2','mismatch_3','allow_mismatches','Description']]
        print('ok.')

    def show_mismatch(self):
        """print mismatch table"""
        return self.mismatch

    def prepare_base_mask(self, total_bases, requested_bases, is_index):
        """format base mask"""
        if requested_bases == 0:
            return 'n*'
        if is_index:# is index read
            if total_bases == requested_bases:
                return 'i*'
            else:
                return 'i{}n{}'.format(requested_bases,total_bases-requested_bases)
        else:# sample read
            if total_bases <= requested_bases:
                return 'y*'
            else:
                return 'y{}n{}'.format(requested_bases,total_bases-requested_bases)

    def infer_base_mask(self, index_len, index2_len, request_seq_type):
        """infer base mask based on sequenced read length and info from samplesheet"""
        # check if index length is longer than the sequenced index bases?
        if index_len > self.index_1_len or index2_len > self.index_2_len:
            print('failed.\nIndex length in samplesheet ({},{}) is greater than the sequenced length ({},{}).'.format(index_len,index2_len,self.index_1_len,self.index_2_len))
            sys.exit(12)
        # get read bases from the 'Description' requested sequencing type
        # infer mask
        mask = []
        if request_seq_type:# non-empty string
            tpat = search('([^\d]+)(\d+)', request_seq_type)
            if not tpat:# unknown pattern
                print('failed.\nFailed to get requested sequencing read bases')
                sys.exit(13)
            sp,rlen = tpat.groups()
            # add read 1
            if self.read_1_len > 0:
                mask.append(self.prepare_base_mask(self.read_1_len, int(rlen), False))
            # add index 1
            if self.index_1_len > 0:
                mask.append(self.prepare_base_mask(self.index_1_len, index_len, True))
            # add index 2
            if self.index_2_len > 0:
                mask.append(self.prepare_base_mask(self.index_2_len, index2_len, True))
            # add read 2
            if self.read_2_len > 0:
                mask.append(self.prepare_base_mask(self.read_2_len, int(rlen), False))
        else:# empty string
            # add read 1
            if self.read_1_len > 0:
                mask.append(self.prepare_base_mask(self.read_1_len, self.read_1_len, False))
            # add index 1
            if self.index_1_len > 0:
                mask.append(self.prepare_base_mask(self.index_1_len, index_len, True))
            # add index 2
            if self.index_2_len > 0:
                mask.append(self.prepare_base_mask(self.index_2_len, index2_len, True))
            # add read 2
            if self.read_2_len > 0:
                mask.append(self.prepare_base_mask(self.read_2_len, self.read_2_len, False))
        return ','.join(mask)

    def write_to_files(self, df):
        """write samplesheet/shell script (function to apply to DataFrame)"""
        # extract info
        index_len = df['index_len'].iloc[0]
        index2_len = df['index2_len'].iloc[0]
        request_seq_type = df['Description'].iloc[0]
        sn = df['sn'].iloc[0]
        mask = self.infer_base_mask(index_len, index2_len, request_seq_type)
        allow_mismatch = df['allow_mismatches'].iloc[0]
        # write to file
        with open('{}/samplesheet_grp{}.csv'.format(self.samplesheet_folder,sn+1),'w') as fsp, open('{}/run_bcl2fq_grp{}.sh'.format(self.script_folder,sn+1), 'w') as frs:
            # samplesheet
            fsp.write(''.join(self.header))
            df.to_csv(fsp, sep=',', columns=self.orig_columns, index=False)
            # shell script
            frs.write('nohup {}  --runfolder-dir {}/  --output-dir {}/Unaligned_{}/'.format(self.bcl2fastq, self.run_folder, self.run_folder, sn+1))
            frs.write('  --sample-sheet {}/samplesheet_grp{}.csv'.format(self.samplesheet_folder, sn+1))
            frs.write('  --use-bases-mask {}'.format(mask))
            if not allow_mismatch:
                frs.write('  --barcode-mismatches 0')
            if self.no_lane_split:
                frs.write('  --no-lane-splitting')
            frs.write('  >run_bcl2fq_grp{}.log\n'.format(sn+1))

    def generate_scripts(self):
        """generate sub-samplesheet and shell script"""
        print('Generating samplesheet and shell script...', end='')
        # add a column indicating the serial number of groups
        self.table = self.table.merge(self.table[['index_len','index2_len','allow_mismatches','Description']].drop_duplicates().reset_index(drop=True).reset_index().rename(columns={'index':'sn'}), how='left', on=['index_len','index2_len','allow_mismatches','Description'])
        # check existence of output files if not allow overwrite
        if not self.overwrite:
            for sn in self.table['sn'].unique():
                samplesheet_file = '{}/samplesheet_grp{}.csv'.format(self.samplesheet_folder,sn+1)
                script_file = '{}/run_bcl2fq_grp{}.sh'.format(self.script_folder,sn+1)
                if exists(samplesheet_file):
                    print('failed.\nsamplesheet file already exists, use -f to allow overwriting.\n{}'.format(samplesheet_file))
                    sys.exit(14)
                if exists(script_file):
                    print('failed.\nscript file already exists, use -f to allow overwriting.\n{}'.format(script_file))
                    sys.exit(15)
        # for each group, create samplesheet and shell script
        self.table.groupby(by=['index_len','index2_len','allow_mismatches','Description']).apply(self.write_to_files)
        print('ok.')

# functions
def get_arguments():
    """fetch commandline arguments."""
    parser = ArgumentParser(description="""Given a samplesheet, check its format and automatically generate shell script for demux.""", 
                            prog='auto_bcl2fq_script.py')
    parser.add_argument("-v", "--version", action="version", version='%(prog)s v2.0')
    parser.add_argument("-r", "--runfolder", nargs="?", required=False, default=getcwd(), help="sequencing run folder", metavar="run_folder", dest="runfolder")
    parser.add_argument("-i", "--samplesheet", nargs="?", required=True, help="samplesheet file", metavar="samplesheet_file", dest="infile")
    parser.add_argument("-s", "--skiplanes", nargs="*", required=False, default=[], type=int, choices=[1,2,3,4,5,6,7,8], help="skip lanes (e.g. -s 2 3 4)", metavar="lanes_to_skip", dest="skip")
    parser.add_argument("-b", "--bcl2fastq", nargs="?", required=False, default="/usr/local/bin/bcl2fastq", help="location of bcl2fastq program", metavar="bcl2fastq", dest="bcl2fastq")
    parser.add_argument("-o", "--scriptdir", nargs="?", required=False, default=None, help="output folder to write shell script(s) (default: run_folder)", metavar="script_folder", dest="scriptdir")
    parser.add_argument("-e", "--samplesheetdir", nargs="?", required=False, default=None, help="output folder to write samplesheet(s) (default: basecall_folder)", metavar="samplesheet_folder", dest="samplesheetdir")
    parser.add_argument("-f", "--force", action="store_true", required=False, default=False, help="whether or not to overwrite output script if existing", dest="force")
    parser.add_argument("-n", "--no-lane-splitting", action="store_true", required=False, default=False, help="whether or not to add --no-lane-splitting option to shell script", dest="nosplit")
    parser.add_argument("-p", "--platform", nargs="?", required=False, default="novaseq", choices=["hiseq","nextseq500","nextseq2000","novaseq"], help="sequencing platform", dest="platform")
    return parser.parse_args()

def main():
    """call me to get started!"""
    args = get_arguments()
    d = MyDemux(args)
    d.get_read_length()
    #d.show_read_length()
    d.load_samplesheet()
    #print(d.table.head())
    d.examine_samplesheet()
    d.calculate_mismatch()
    #print(d.show_mismatch())
    d.generate_scripts()

# main
if __name__ == '__main__':
    main()

