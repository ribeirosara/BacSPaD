# Quick start guide


Welcome to BacSPaD, your comprehensive resource for bacterial strains pathogenicity data. Here’s a brief overview to help you get started:

**Data Section**

Perform keyword queries across all fields or specific fields. For instance, if you want to focus your analysis on certain species, genus, families, order, class or phylum select  the corresponding metadata fields. You may also be eventually interested in selecting certain diseases using “disease category”, “disease subcategory”, and “isolation source” fields.

Use BV-BRC genome IDs to retrieve corresponding genomes and data files (e.g., proteomes, protein families) from the BV-BRC FTP website (https://www.bv-brc.org/docs/quick_references/ftp.html) as there demonstrated:

“Once you have copied the list of genome ids you are interested in a separate file “genome_list”, you can use the following one line shell script to read the list of genome ids from your file and download corresponding .fna files from the PATRIC FTP site. If you are interested in other file type, say .PATRIC.faa or .PATRIC.features.tab, simply replace .fna with that extension.


<code>
for i in \` cat genome_list \`; do wget -qN "ftp://ftp.bvbrc.org/genomes/\$i/\$i.fna";
done
</code>

**Dashboard**

Interactive Visualizations: Explore top species and families, taxa distribution, associated disease categories and isolation sources, and a location map of all strains.

**Molecular Biology**

Genome Data: Visualize distributions of the number of plasmids, contigs, genome lengths, GC content (%), and number of coding sequences for both HP and NHP genomes.

**Virulence Factors**

Access information on virulence factors for the most prevalent clinical species, including gene name, gene frequency in HP and NHP strains, the associated BV-BRC genome ids and species names; as well as the number of associated strains, species, genus, and families.

**About**

Learn about BacSPaD’s utility for microbiology research and find contact information.

**Citation**

If you use BacSPaD, please cite it as follows:

Ribeiro, S., Chaumet, G., Alves, K., Nourikyan, J., Shi, L., Lavergne,J.-P., Mijakovic, I., de Bernard, S., & Buffat, L. (2024). BacSPaD: A robust bacterial strains’ pathogenicity resource based on integrated and curated genomic metadata. Preprints, 202407.0837.v1. 
https://doi.org/10.20944/preprints202407.0837.v1


For further assistance, contact us at sara.ribeiro@altrabio.com.

