'''
Date: 2021-02-04 15:53:53
LastEditors: cany
LastEditTime: 2021-02-22 16:52:50
'''

from anarci import number
import pandas as pd
import numpy as np
import os
import json
import re
from pandas.io.json import json_normalize
import pysnooper

__author__ = 'yincan'
__contact__ = 'yin_can@wuxibiologics.com'


class divideScheme:
    '''
    Almost base on IMGT numbering scheme.
    The boundary of CDR2 broadens by 1 to the left and right. 
    '''
    def __init__(self):
        self.getCDR()

    @staticmethod
    def convertRange(start, end):
        extent = [i for i in range(start, end + 1)]
        str_ext = [st for st in map(str, extent)]
        return str_ext

    def getCDR(self):
        self.cdr1 = self.convertRange(27, 38)
        #The boundary of CDR2 broadens by 1 to the left and right.
        self.cdr2 = self.convertRange(55, 66)
        self.cdr3 = self.convertRange(105, 117)


class Fasta2Library:
    def __init__(self, fasta, scheme='imgt'):
        self.fasta = fasta
        self.scheme = scheme
        self.out = self.fasta.rsplit(".")[0]

    ######################### Basic Function ########################
    def seqS(self):
        '''
        Merge multi-line fasta file to one-line fasta sequence.
        All the fastas consist of seqlist.
        '''
        seqlist = []
        with open(self.fasta, 'r') as f:
            seq = ''
            for line in f:
                if not line.startswith('>'):
                    seq += line.replace('\n', '').strip()
                else:
                    if seq != '':
                        seqlist.append(seq)
                    seq = ''
        return seqlist

    def outNumbering(self):
        '''
        Process on number output of numbering
        '''
        dict_numbering = {}
        for tmp in self.numbering:
            pos = str(tmp[0][0])
            AA = tmp[1]
            dict_numbering[pos] = AA

        return dict_numbering

    @staticmethod
    def natural_sort_key(s):
        _nsre = re.compile('([0-9]+)')
        return [
            int(text) if text.isdigit() else text.lower()
            for text in re.split(_nsre, s)
        ]

    @staticmethod
    def _getVarRg(varglist, mycolumns):
        '''
        Extract columns of varible region, get pos from mydata.
        Mycolumns conclude 111A ,111B such as the non-numeric class.
        '''
        varpos = [
            pos for pos in mycolumns
            if re.search('([0-9]+)', pos).group() in varglist
        ]

        return varpos

    def extractRg2df(self, df, Rg):
        '''
        经过处理或这原始data(df)，提取需要的区域的pos，再从df中提取这部分dataframe(RgDf)
        '''
        self.mycolumns = self.df.columns.tolist(
        )  # mycolumns conclude 111A ,111B such as the non-numeric class

        Pos = self._getVarRg(Rg, self.mycolumns)
        RgDf = df[Pos].round(decimals=3)

        return RgDf

    # @pysnooper.snoop()
    def divideReg(self, df):
        self.CDR = divideScheme()
        self.cdr1df = self.extractRg2df(df, self.CDR.cdr1)
        self.cdr2df = self.extractRg2df(df, self.CDR.cdr2)
        self.cdr3df = self.extractRg2df(df, self.CDR.cdr3)

    #########################    Exec    ########################
    def fasta2df(self):
        '''
        Numbering sequence by ANARCI's function "number".
        Seqs in fasta file convert to dataframe.
        '''
        seqlist = self.seqS()
        self.data = []
        for self.seq in seqlist:
            try:
                self.numbering, chain_type = number(self.seq, self.scheme)
            except AssertionError:
                print(self.seq)
            pos_aaDic = self.outNumbering()
            self.data.append(pos_aaDic)
        df = json_normalize(self.data)
        self.df = df.reindex(sorted(df.columns, key=self.natural_sort_key),
                             axis=1).replace("-", np.nan)

    def df2Library(self):
        '''
        Calculate distribution.
        Attention!!!!! 
        The order of row and column!!!!        
        '''
        self.fasta2df()
        self.div = self.df.apply(lambda x: pd.value_counts(x, normalize=True),
                                 axis=0).fillna(0)
        '''
        Dataframe is divide into cdrs.
        '''
        self.divideReg(self.div)
        self.cdr1df.T.to_csv(f'{self.out}.cdr1.csv', index_label='pos')
        self.cdr2df.T.to_csv(f'{self.out}.cdr2.csv', index_label='pos')
        self.cdr3df.T.to_csv(f'{self.out}.cdr3.csv', index_label='pos')

    def diffLength(self):
        self.fasta2df()
        self.divideReg(self.df)
        for num, cdr in enumerate([self.cdr1df, self.cdr2df, self.cdr3df]):
            name = 'cdr' + str(num + 1)
            cdr['len'] = cdr.count(axis=1)
            vdf = cdr.len.value_counts(normalize=True)
            lendf = pd.DataFrame({'length': vdf.index, 'percentage': vdf.values})
            lendf.to_csv(f'LengthDistribution.{self.out}.{name}.csv',index=None)

            for length in vdf.index:
                pop_df = cdr[cdr['len'] == length]

                pop_df.drop(['len'], axis=1, inplace=True)
                mydf = pop_df.dropna(axis=1, how='all')
                ww_df = mydf.apply(
                    lambda x: pd.value_counts(x, normalize=True), axis=0)

                ff_df = ww_df.fillna(0).T

                ff_df.to_csv(f'{self.out}.{name}.len{length}.csv',
                             index_label='pos')


class OutCdr(Fasta2Library):
    def __init__(self, fasta, scheme='imgt'):
        Fasta2Library.__init__(self, fasta, scheme='imgt')
        self.cvtFas()
        self.outF()

    def cvtEseq(self, seq):
        '''
        Convert sequence to dataframe for each seq.
        '''
        self.numbering, chain_type = number(seq, self.scheme)
        pos_aaDic = self.outNumbering()
        seqdf = pd.DataFrame([pos_aaDic])
        self.seqdf = seqdf.reindex(sorted(seqdf.columns,
                                          key=self.natural_sort_key),
                                   axis=1)

    def cvtFas(self):
        seqlist = self.seqS()
        self.data = []
        for self.seq in seqlist:
            try:
                self.numbering, chain_type = number(self.seq, self.scheme)
            except AssertionError:
                print(self.seq)
            pos_aaDic = self.outNumbering()
            self.data.append(pos_aaDic)
        df = json_normalize(self.data)
        self.df = df.reindex(sorted(df.columns, key=self.natural_sort_key),
                             axis=1).replace(np.nan, "-")

    def outF(self):
        self.divideReg(self.df)
        self.cdr1df.apply(lambda x: ''.join(x),
                          axis=0).drop_duplicates().to_csv(
                              f'{self.out}.cdr1.csv', index=None)
        self.cdr2df.apply(lambda x: ''.join(x),
                          axis=0).drop_duplicates().to_csv(
                              f'{self.out}.cdr2.csv', index=None)
        self.cdr3df.apply(lambda x: ''.join(x),
                          axis=0).drop_duplicates().to_csv(
                              f'{self.out}.cdr3.csv', index=None)


if __name__ == "__main__":
    import argparse
    import time
    t1 = time.time()
    parser = argparse.ArgumentParser(
        prog='Fasta2Library',
        description='\t\tConvert fasta to a.a. distribution',
        epilog=('contact: {} <{}>'.format(__author__, __contact__)),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        '-i',
        '--fasta',
        help='Fasta file is necessary, which includes the target sequences',
        required=True)
    parser.add_argument('-n',
                        '--scheme',
                        choices=('imgt', 'kabat'),
                        default='imgt',
                        help='Antibody numbering scheme',
                        required=False)
    parser.add_argument('-d',
                        '--divide',
                        action='store_true',
                        help='Output divided sequences of cdr regions.',
                        required=False)
    parser.add_argument(
        '-l',
        '--length',
        action='store_true',
        help='Output a.a. distribution with different length of cdr regions.',
        required=False)

    args = vars(parser.parse_args())
    fasta = args.get('fasta')
    scheme = args.get('scheme')
    if args.get('divide'):
        OutCdr(fasta=fasta, scheme=scheme)
    elif args.get('length'):
        Fasta2Library(fasta=fasta, scheme=scheme).diffLength()
    else:
        Fasta2Library(fasta=fasta, scheme=scheme).df2Library()

    t2 = time.time()
    print(t2 - t1)