select * from snp135 where name='rs10019843';

select name, chrom, count(*) as mappings
from snp135 s
group by name, chrom
having mappings > 1
limit 10;

select * from snp135 limit 100;
