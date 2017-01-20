"""
Partition data sets by host, annotate with collection dates
"""
from datetime import datetime
from seqUtils import convert_fasta
from csv import DictReader
from glob import glob

def process(fasta_file, csvdict):
    handle = open(fasta_file, 'rU')
    fasta = dict(convert_fasta(handle))
    handle.close()

    handle = open(fasta_file.replace('.fa', '.gb.csv'), 'rU')
    rows = DictReader(handle)
    for row in rows:
        # header,accno,patid,coldate,tissue,moltype,isolate,isosource
        h = row['header']
        try:
            coldate = datetime.strptime(row['coldate'], '%d-%b-%Y')
        except:
            print row
            raise

        days = (coldate-datetime(2000,1,1)).days
        patid = row['patid']
        if patid not in csvdict:
            csvdict.update({patid: {}})

        csvdict[patid].update({h: {
            'accno': row['accno'],
            'isolate': row['isolate'],
            'moltype': row['moltype'].split()[-1],
            'days': days,
            'tissue': row['tissue'],
            'source': row['isosource'],
            'seq': fasta[h]
        }})


def write_output(stem, d):
    for patid, data in d.iteritems():
        outfile = open('%s.%s.fa' % (stem, patid), 'w')
        for h, annot in data.iteritems():
            try:
                outfile.write('>%s_%s_%s_%s_%s_%s_%d\n%s\n' % (
                    annot['accno'],
                    patid,
                    annot['isolate'],
                    annot['moltype'],
                    annot['tissue'],
                    annot['source'],
                    annot['days'],
                    annot['seq']
                ))
            except:
                print annot
                raise
        outfile.close()


def main():
    files = filter(lambda x: not x.endswith('.annot.fa'),
                   glob('../data/25712966/pol-*.fa'))
    csvdict = {}
    for file in files:
        process(file, csvdict)
    write_output('../data/25712966/pol.annot', csvdict)

if __name__ == '__main__':
    main()
