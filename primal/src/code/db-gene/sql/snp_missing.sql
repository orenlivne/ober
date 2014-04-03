#----------------------------------------------------------------
# Find Hutt Affy SNPs that map to multiple chromosomes.
# Run on local UCSC+SNP database
#
# Author: Oren Livne
# Date: 15-NOV-2012
#----------------------------------------------------------------

# Find Hutt SNP RS#s that have multiple chromosome matches in dbSNP135; save results to file
# Warning: this query will be slow if over 60,000 records (maybe 65,536?) in MySQL on my
# home desktop. That's why we break it into parts. This query is called within a bash wrapper
# script loop.

# Load the list of Hutterites SNPs (exported from plink) into a table
drop table if exists snp_hutt;
create temporary table snp_hutt (snp varchar(20));
delete from snp_hutt;
    
# Must be a literal string
load data infile 'c:/ober/data/common/db/snp_missing.nof' 
into table snp_hutt lines terminated by '\r\n';

create index snp_hutt_snp on snp_hutt (snp);
        
select m.snp, count(distinct s.chrom) as chrom_mappings 
from 
    snp_hutt m 
    left join snp135 s on s.name = m.snp 
group by m.snp 
having chrom_mappings > 1
into outfile 'c:/ober/data/common/db/snp_missing.out';
