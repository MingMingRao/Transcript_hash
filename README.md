# Transcript_hash
This script（transcript_hash.pl） converts the known HGVS information to the coordinates on the hg19. You just need to provide the HGVS you want to convert and GRCh37_latest_genomic.gff（from ncbi）.The format can be as follows：
1、A file with suffix GTF or xls
line1   nonsynonymous SNV       MTOR:NM_004958:exon40:c.C5664A:p.F1888L,        1 11189845 11189845 G T
line2   nonsynonymous SNV       MTOR:NM_004958:exon30:c.G4448T:p.C1483F,        1 11217230 11217230 C A
line4   nonsynonymous SNV       MPL:NM_005373:exon10:c.T1543G:p.W515G,  1 43815008 43815008 T G
line5   nonsynonymous SNV       MPL:NM_005373:exon10:c.G1544C:p.W515S,  1 43815009 43815009 G C
line6   nonsynonymous SNV       NRAS:NM_002524:exon3:c.C176A:p.A59D,    1 115256535 115256535 G T
2、A file with suffix txt
NM_004972       c.A2047G
NM_005157       c.A1187G
NM_000222       c.T2466G
NM_007294       c.A3113G
3、Direct information like this format: NM_004972:c.A2047G is entered on the first parameter of the script

Usage method：
   perl transcript_hash.pl input GRCh37_latest_genomic.gff output（The input file format can be more than three
   The output file can be you change）
