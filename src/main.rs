use clap::{App, load_yaml};
use rust_htslib::bcf::record::GenotypeAllele;
use rust_htslib::bcf::{Reader, Read};
use std::collections::HashMap;
use std::collections::hash_map::Entry;
use rayon::prelude::*;

fn main() {
    let yaml = load_yaml!("cli.yml");
    let matches = App::from_yaml(yaml).get_matches();

    let vcf = matches.value_of("vcf").unwrap();

    let mut bcf = Reader::from_path(vcf).expect("Error opening file.");

    let contigs = vec!["12"];

    let mut cont_hm_postn_vec_geno_vec: HashMap<String, (Vec<usize>, Vec<Vec<Vec<GenotypeAllele>>>)> = HashMap::new();

    for record_result in bcf.records() {
        let record = record_result.as_ref().expect("Failed to read record!");
        let desc = record.desc();
        let cont_pos: Vec<&str> = desc.split(":").collect();
        if contigs.iter().any(|x| *x == cont_pos[0]) {
            let sample_count = record.sample_count().to_string().parse::<usize>().unwrap();
            let cont = cont_pos[0].to_string();
            let postn = cont_pos[1].to_string().parse::<usize>().unwrap();
            let genos = record.genotypes().expect("Failed to read genotypes");
            let simple_genos: Vec<Vec<GenotypeAllele>> = (0..sample_count)
                .map(|sample_index| {
                    let sample_geno: Vec<GenotypeAllele> = genos.get(sample_index).iter().map(|allele| *allele).collect();
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

    let cont_hm_post_maf_vec: HashMap<String, Vec<(usize, f64)>> = cont_hm_postn_vec_geno_vec
        .iter()
        .map(|(cont, (postn_vec, geno_vec))| {
            let postn_maf: Vec<(usize, f64)> = postn_vec
                .par_iter()
                .zip(geno_vec)
                .map(|(postn, genos)| {
                    let geno_counts = genos
                        .iter()
                        .fold((0, 0, 0), |mut acc, geno| {
                            for allele in geno.iter() {
                                match allele.index() {
                                    Some(0) => {
                                        acc.0 += 1usize;
                                    },
                                    Some(_) => {
                                        acc.1 += 1usize;
                                    },
                                    None => {
                                        acc.2 += 1;
                                    },
                                }
                            }
                            acc
                        });

                    let num_alleles = (postn_vec.len() * 2) as f64;
                    let num_uncalled = geno_counts.2 as f64;
                    let num_called_alleles = num_alleles - num_uncalled;
                    let maf = geno_counts.0.min(geno_counts.1) as f64 / num_called_alleles;
                    (postn.clone(), maf.clone())
                })
                .filter(|(_, maf)| *maf > 0.05)
                .collect();
            println!("{:?}", postn_maf);
            (cont.clone(), postn_maf.clone())
        })
        .collect();

    


    // Iter through MAFs (use enumerate) and calc r2 for up to N SNPs upstream and downstream of each SNP pair, calc SNP "Sparsity"
    // Take least sparse SNP and perform pruning using HDBSCAN clustering and distance threshold
    // Recalculate Sparsities of neighbouring SNPs
    // repeat until provided fraction remain
}
