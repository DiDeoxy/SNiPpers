use clap::{App, load_yaml};
use rust_htslib::bcf::{Reader, Read};
use std::collections::HashMap;

fn main() {
    let yaml = load_yaml!("cli.yml");
    let matches = App::from_yaml(yaml).get_matches();

    let vcf = matches.value_of("vcf").unwrap();

    let mut bcf = Reader::from_path(vcf).expect("Error opening file.");

    let contigs = vec!["1"];

    let cont_pos_geno_vec: Vec<(String, usize, Vec<Vec<i32>>)> = bcf
        .records()
        .filter(|record_result| {
            let record = record_result.as_ref().expect("Failed to read record");
            let desc = record.desc();
            let cont_pos: Vec<&str> = desc.split(":").collect();
            if contigs.iter().any(|x| x == &cont_pos[0]) {
                return true
            } else {
                return false
            }
        })
        .map(|record_result| {
            let record = record_result.expect("Failed to read record");
            let sample_count = record.sample_count().to_string().parse::<usize>().unwrap();
            let desc = record.desc();
            let cont_pos: Vec<&str> = desc.split(":").collect();
            let pos = cont_pos[1].to_string().parse::<usize>().unwrap();
            let genos = record.genotypes().expect("failed to read genotypes");
            let simple_genos: Vec<Vec<i32>> = (0..sample_count)
                .map(|sample_index| {
                    let sample_geno: Vec<i32> = genos.get(sample_index).iter().map(|allele| i32::from(*allele)).collect();
                    sample_geno
                })
                .collect();
            return (cont_pos[0].to_string(), pos, simple_genos)
        })
        .collect();

    // for record_result in bcf.records() {
    //     let record = record_result.expect("Failed to read record");
    //     let desc = record.desc();
    //     let cont_pos: Vec<&str> = desc.split(":").collect();
    //     // let genos = record.genotypes().expect("failed to read genotypes");

    //     println!("{:?}", std::str::from_utf8(&record.id()));
    //     println!("{:?}", cont_pos);

    //     break
    // }
}
