use clap::{App, load_yaml};
use rust_htslib::bcf::{Reader, Read};
use std::collections::HashMap;
use std::collections::hash_map::Entry;

fn main() {
    let yaml = load_yaml!("cli.yml");
    let matches = App::from_yaml(yaml).get_matches();

    let vcf = matches.value_of("vcf").unwrap();

    let mut bcf = Reader::from_path(vcf).expect("Error opening file.");

    let contigs = vec!["1"];

    let mut cont_hm_postn_vec_geno_vec: HashMap<String, (Vec<usize>, Vec<Vec<Vec<i32>>>)> = HashMap::new();

    for record_result in bcf.records() {
        let record = record_result.as_ref().expect("Failed to read record!");
        let desc = record.desc();
        let cont_pos: Vec<&str> = desc.split(":").collect();
        if contigs.iter().any(|x| *x == cont_pos[0]) {
            let sample_count = record.sample_count().to_string().parse::<usize>().unwrap();
            let postn = cont_pos[1].to_string().parse::<usize>().unwrap();
            let cont = cont_pos[0].to_string();
            let genos = record.genotypes().expect("Failed to read genotypes");
            let simple_genos: Vec<Vec<i32>> = (0..sample_count)
                .map(|sample_index| {
                    let sample_geno: Vec<i32> = genos.get(sample_index).iter().map(|allele| i32::from(*allele)).collect();
                    sample_geno
                })
                .collect();

            let (postn_vec, geno_vec) = match cont_hm_postn_vec_geno_vec.entry(cont) {
                Entry::Occupied(entry) => entry.into_mut(),
                Entry::Vacant(entry) => entry.insert((Vec::new(), Vec::new()))
            };

            postn_vec.push(postn);
            geno_vec.push(simple_genos);
        }
    }



    // Iter through snps, calc MAF, filter out MAFs below threshhold
    // Iter through MAFs (use enumerate) and calc r2 for up to N SNPs upstream and downstream of each SNP pair, calc SNP "Sparsity"
    // Take least sparse SNP and perform pruning using HDBSCAN clustering and distance threshold
    // Recalculate Sparsities of neighbouring SNPs
    // repeat until provided fraction remain
}
