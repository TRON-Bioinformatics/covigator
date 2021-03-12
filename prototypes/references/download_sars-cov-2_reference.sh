wget ftp://ftp.ensemblgenomes.org/pub/viruses/fasta/sars_cov_2/cdna/*
wget ftp://ftp.ensemblgenomes.org/pub/viruses/fasta/sars_cov_2/dna/*
wget ftp://ftp.ensemblgenomes.org/pub/viruses/fasta/sars_cov_2/pep/*
wget ftp://ftp.ensemblgenomes.org/pub/viruses/gff3/sars_cov_2/*
wget ftp://ftp.ensemblgenomes.org/pub/viruses/json/sars_cov_2/*

gzip -d *.fa.gz
gzip -d *.gff3.gz