#----------------------------------------------------------------
# Find Affy SNPs that map to multiple chromosomes in dbSNP 135.
# Run against our local UCSC+SNP database.
#
# Author: Oren Livne
# Date: 16-NOV-2012
#----------------------------------------------------------------

# Find Hutt SNP RS#s that have multiple chromosome matches in dbSNP135; save results to file
# Warning: this query will be slow if over 60,000 records (maybe 65,536?) in MySQL on my
# home desktop. That's why we break it into parts. This query is called within a bash wrapper
# script loop.

# Load the list of Hutterites SNPs (exported from plink) into a table
drop table if exists snp_list;
create temporary table snp_list (snp varchar(20));
delete from snp_list;
    
# Must be a literal string
load data infile '${INPUT_FILE}' 
into table snp_list lines terminated by '\r\n';

create index snp_list_snp on snp_list (snp);

select m.snp, count(distinct s.chrom) as chrom_mappings 
from 
    snp_list m 
    left join snp135 s on s.name = m.snp 
group by m.snp 
having chrom_mappings > 1
into outfile '${OUTPUT_FILE}';
