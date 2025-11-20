Before starting to write code for the request below, formulate any questions that you still have for me in order to implement the request correctly.

add a new subcommand, `inquiSTR benchmark`, that takes in 1) a inquiSTR call output file and 2) a VCF file with tandem repeat genotypes and 3) a mode, which is for now either MAX (default) or MIN and 4) a --plot argument determining an output file.

The VCF file and inquiSTR call output are assumed to contain genotypes for the same repeats, and those that are only found in one of the files will be dropped. For now we will assume that variants can be matched based on the chromosome and START/FIRST coordinate of the inquiSTR call file, and the chromosome and POS field of the VCF file. The END is to be ignored.
From the VCF file, the length of the ALT sequence(s) determines the repeat genotype.  The output of the subcommand is a correlation scatter plot of inquiSTR lengths vs. VCF lengths, picking either the longest allele per locus (MAX; the default) or the shortest allele per locus (MIN). The plot is generated using plotly.

I would suggest first iterating over the inquiSTR call file, and then over the VCF taking from the VCF those chrom:pos records that match with an allele in the inquiSTR call file.

Also print the Pearson correlation coefficient, the number of variants in the inquiSTR call file found in the VCF file (and vice versa) and the number of variants found in both.