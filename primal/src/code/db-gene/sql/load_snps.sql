# Load the list SNps into a temporary table
drop table if exists snp_temp;
create temporary table snp_temp (snp varchar(20));
delete from snp_temp;

# Get SNP metadata from the UCSC table and store in our table
delete from snp;
select snp_temp.snp, count(distinct s.chrom) as chrom_mappings 
from 
    snp_temp
    left join snp135 s on s.name = m.snp 
group by m.snp 
