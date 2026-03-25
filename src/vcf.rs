use std::io::{self, Write};

struct SampleFields {
    gt: &'static str,
    rb: [String; 2],
    frb: [String; 2],
    svlen: String,
}

fn sample_fields(phase1: Option<f64>, phase2: Option<f64>, interval_length: f64) -> SampleFields {
    let fmt = |v: f64| format!("{:.0}", v);
    let missing = ".".to_string();

    match (phase1, phase2) {
        (None, None) => SampleFields {
            gt: "./.",
            rb: [missing.clone(), missing.clone()],
            frb: [missing.clone(), missing.clone()],
            svlen: missing,
        },
        (None, Some(p2)) => SampleFields {
            gt: "./1",
            rb: [missing.clone(), fmt(p2)],
            frb: [missing.clone(), fmt(p2 + interval_length)],
            svlen: fmt(p2),
        },
        (Some(p1), None) => SampleFields {
            gt: "1/.",
            rb: [fmt(p1), missing.clone()],
            frb: [fmt(p1 + interval_length), missing],
            svlen: fmt(p1),
        },
        (Some(p1), Some(p2)) => {
            let gt = if (p1 - p2).abs() < 0.5 { "1|1" } else { "1|2" };
            SampleFields {
                gt,
                rb: [fmt(p1), fmt(p2)],
                frb: [fmt(p1 + interval_length), fmt(p2 + interval_length)],
                svlen: fmt(p1.max(p2)),
            }
        }
    }
}

/// Write VCF to stdout by reading from a TSV file (avoids keeping all genotypes in memory)
pub fn write_vcf_to_stdout(
    tsv_path: &std::path::Path,
    sample: &str,
    reference: &Option<String>,
    bam: &str,
) {
    use std::fs::File;
    use std::io::{BufRead, BufReader};

    let stdout = io::stdout();
    let mut writer = io::BufWriter::new(stdout);

    // Write VCF header
    writeln!(writer, "##fileformat=VCFv4.3").expect("Failed writing VCF header");
    writeln!(writer, "##source=inquiSTR").expect("Failed writing VCF header");

    // Add reference if provided
    if let Some(ref_path) = reference {
        writeln!(writer, "##reference={}", ref_path).expect("Failed writing VCF header");
    }

    // Add contig lines from BAM/CRAM header
    if let Ok(chrom_lengths) =
        crate::bam_utils::get_chrom_lengths_from_bam_header(bam.to_string(), reference)
    {
        let mut contigs: Vec<_> = chrom_lengths.into_iter().collect();
        contigs.sort_by(|a, b| a.0.cmp(&b.0));
        for (name, length) in contigs {
            writeln!(writer, "##contig=<ID={},length={}>", name, length)
                .expect("Failed writing VCF header");
        }
    }

    // Add INFO fields
    writeln!(
        writer,
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant\">"
    )
    .expect("Failed writing VCF header");
    writeln!(
        writer,
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">"
    )
    .expect("Failed writing VCF header");
    writeln!(writer, "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">").expect("Failed writing VCF header");

    // Add FORMAT fields
    writeln!(writer, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
        .expect("Failed writing VCF header");
    writeln!(
        writer,
        "##FORMAT=<ID=RB,Number=2,Type=Integer,Description=\"Repeat length of the two alleles in bases relative to reference\">"
    )
    .expect("Failed writing VCF header");
    writeln!(
        writer,
        "##FORMAT=<ID=FRB,Number=2,Type=Integer,Description=\"Full repeat length of the two alleles in bases\">"
    )
    .expect("Failed writing VCF header");

    // Add ALT definition for STR
    writeln!(writer, "##ALT=<ID=STR,Description=\"Short Tandem Repeat\">")
        .expect("Failed writing VCF header");

    // Write column headers
    writeln!(writer, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}", sample)
        .expect("Failed writing VCF header");

    // Read TSV and convert to VCF records
    let tsv_file = File::open(tsv_path).expect("Failed to open temp TSV file");
    let reader = BufReader::new(tsv_file);

    let mut idx = 0;
    for line in reader.lines() {
        let line = line.expect("Failed reading TSV line");

        // Skip header line
        if line.starts_with("chromosome") || line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != crate::filetype::STR_MIN_COLUMNS {
            eprintln!(
                "ERROR: Malformed line in TSV (expected {} columns, got {})",
                crate::filetype::STR_MIN_COLUMNS,
                fields.len()
            );
            std::process::exit(1);
        }

        let chrom = fields[0];
        let start: u32 = fields[1].parse().expect("Invalid start coordinate in TSV");
        let end: u32 = fields[2].parse().expect("Invalid end coordinate in TSV");
        let _info = fields[3]; // Currently unused, but available if needed
        let phase1_str = fields[crate::filetype::STR_FIXED_COLUMNS];
        let phase2_str = fields[crate::filetype::STR_FIXED_COLUMNS + 1];

        idx += 1;
        let pos = start + 1; // VCF is 1-based
        let id = format!("STR_{}", idx);

        // Parse phase values (handle "nan" or numeric values)
        let phase1: Option<f64> = phase1_str.parse().ok();
        let phase2: Option<f64> = phase2_str.parse().ok();

        let s = sample_fields(phase1, phase2, (end - start) as f64);

        let info = format!("END={};SVTYPE=STR;SVLEN={}", end, s.svlen);
        let sample_data = format!("{}:{},{}:{},{}", s.gt, s.rb[0], s.rb[1], s.frb[0], s.frb[1]);

        writeln!(
            writer,
            "{}\t{}\t{}\tN\t<STR>\t.\t.\t{}\tGT:RB:FRB\t{}",
            chrom, pos, id, info, sample_data
        )
        .expect("Failed writing VCF record");
    }

    writer.flush().expect("Failed to flush VCF to stdout");
}
