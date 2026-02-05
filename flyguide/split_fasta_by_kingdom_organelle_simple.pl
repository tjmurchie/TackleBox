#!/usr/bin/env perl
use strict;
use warnings;

# Split a combined NCBI FASTA (from NCBI-NT_Downloader.pl) into separate files
# by:
#   - Kingdom  (Animal / Plant / Fungi / Bacteria / Archaea / Protist / Protozoa / Unknown)
#   - Type     (Mito / Plastid / NucMark / Other)
#
# Usage:
#   perl split_fasta_by_kingdom_organelle_simple.pl \
#       INPUT.fasta \
#       OUTPREFIX \
#       species_kingdom.tsv
#
# This will produce files like:
#   OUTPREFIX.Animal-Mito.fasta
#   OUTPREFIX.Plant-Plastid.fasta
#   OUTPREFIX.Fungi-NucMark.fasta
#   OUTPREFIX.Unknown-Other.fasta
#
# Assumptions:
#   - FASTA headers are in standard NCBI style, e.g.:
#       >AY425443.1 Americobdella valdiviana cytochrome c oxidase ...
#   - species_kingdom.tsv has:
#       species<TAB>kingdom
#     where "species" is a binomial name (Genus species).

# ----------------- Args -----------------

my ($fasta_in, $out_prefix, $sp_kd_tsv) = @ARGV;
if ( !defined $fasta_in || !defined $out_prefix || !defined $sp_kd_tsv ) {
    die "Usage: $0 INPUT.fasta OUTPREFIX species_kingdom.tsv\n";
}

# ----------------- Load species ? kingdom map -----------------

my %species_to_kd;  # lc(species) -> kingdom
my %genus_to_kd;    # lc(genus)   -> kingdom
my $rows = 0;

open my $tfh, '<', $sp_kd_tsv
    or die "Cannot open species?kingdom TSV [$sp_kd_tsv]: $!\n";

while ( my $line = <$tfh> ) {
    chomp $line;
    next unless $line =~ /\S/;
    my ($sp, $kd) = split /\t/, $line;
    next unless defined $sp && defined $kd;
    $sp =~ s/^\s+//;
    $sp =~ s/\s+$//;
    $kd =~ s/^\s+//;
    $kd =~ s/\s+$//;
    next unless $sp && $kd;

    $rows++;
    my $sp_lc = lc $sp;
    $species_to_kd{$sp_lc} = $kd;

    # Derive genus from species
    my ($genus) = split /\s+/, $sp;
    if ($genus) {
        my $g_lc = lc $genus;
        $genus_to_kd{$g_lc} //= $kd;  # first kingdom wins if multiple
    }
}
close $tfh;

warn "Loaded species_kingdom map: $rows rows, "
   . scalar(keys %species_to_kd) . " species keys, "
   . scalar(keys %genus_to_kd)   . " genus keys\n";

# Show a couple of example lookups if available
foreach my $example (qw(americobdella valdiviana eisenia fetida)) {
    if (exists $species_to_kd{$example}) {
        warn "Example mapping: $example -> $species_to_kd{$example}\n";
    }
}

# ----------------- Helpers -----------------

sub normalize_kingdom_label {
    my ($kd_raw) = @_;
    return 'Unknown' unless defined $kd_raw && $kd_raw ne '';

    my $k = lc $kd_raw;

    return 'Animal'   if $k =~ /animal/;
    return 'Plant'    if $k =~ /plant|plantae/;
    return 'Fungi'    if $k =~ /fungi/;
    return 'Bacteria' if $k =~ /bacteria/;
    return 'Archaea'  if $k =~ /archaea/;
    return 'Protist'  if $k =~ /protist|chromista/;
    return 'Protozoa' if $k =~ /protozoa/;

    return 'Unknown';
}

sub classify_type {
    my ($header_lc) = @_;

    # Priority order: Mito > Plastid > NucMark > Other

    # Mitochondrial markers
    if ( $header_lc =~ /mitochondr/ ||
         $header_lc =~ /\[location=mitochondrion\]/ ) {
        return 'Mito';
    }

    # Plastid / chloroplast
    if ( $header_lc =~ /chloroplast|plastid|plastome/ ||
         $header_lc =~ /\[location=chloroplast\]/ ||
         $header_lc =~ /\[location=plastid\]/ ) {
        return 'Plastid';
    }

    # Nuclear markers: 18S, 28S, SSU, LSU, ITS, 5.8S, histone H3
    if ( $header_lc =~ /\b18s\b/ ||
         $header_lc =~ /\b28s\b/ ||
         $header_lc =~ /small subunit ribosomal|ssu rna/ ||
         $header_lc =~ /large subunit ribosomal|lsu rna/ ||
         $header_lc =~ /\bits1\b|\bits2\b|internal transcribed spacer/ ||
         $header_lc =~ /5\.8s/ ||
         $header_lc =~ /histone h3|\bh3\b/ ) {
        return 'NucMark';
    }

    return 'Other';
}

sub parse_accession {
    my ($header) = @_;
    # Remove leading '>'
    $header =~ s/^>//;

    # Common NCBI format: first token is accession
    my @fields = split /\s+/, $header;
    return $fields[0] if @fields;

    return undef;
}

sub guess_species_name {
    my ($header) = @_;
    # Remove leading '>'
    $header =~ s/^>//;

    my @fields = split /\s+/, $header;
    # skip if too short
    return undef if @fields < 3;

    # accession = fields[0], then Genus = fields[1], species epithet = fields[2]
    my $genus  = $fields[1];
    my $sp_ep  = $fields[2];

    # Simple sanity check: genus starts with capital letter?
    return undef unless $genus =~ /^[A-Z]/;

    return "$genus $sp_ep";
}

# ----------------- Main splitting -----------------

open my $ffh, '<', $fasta_in
    or die "Cannot open FASTA file [$fasta_in]: $!\n";

my %seen_acc;
my %acc_category;
my %fh_for_cat;
my %cat_counts;
my $dup_count        = 0;
my $conflict_count   = 0;
my $current_header   = undef;
my $current_seq      = '';

sub flush_record {
    my ($header, $seq) = @_;
    return unless defined $header && $header ne '';

    my $acc = parse_accession($header);
    return unless $acc;

    # Deduplicate by accession
    if (exists $seen_acc{$acc}) {
        # Check for classification conflict if we wanted to:
        #   (but we skip record anyway)
        $dup_count++;
        return;
    }

    my $header_lc = lc $header;

    # Species / kingdom assignment
    my $species_guess = guess_species_name($header);
    my $kd_raw;
    my $kingdom_label = 'Unknown';

    if ($species_guess) {
        my $sp_lc = lc $species_guess;
        if (exists $species_to_kd{$sp_lc}) {
            $kd_raw = $species_to_kd{$sp_lc};
        } else {
            # genus-only fallback
            my ($genus) = split /\s+/, $species_guess;
            if ($genus) {
                my $g_lc = lc $genus;
                $kd_raw = $genus_to_kd{$g_lc} if exists $genus_to_kd{$g_lc};
            }
        }
    }

    $kingdom_label = normalize_kingdom_label($kd_raw);

    # Genome/marker type
    my $type_label = classify_type($header_lc);

    # If kingdom is still Unknown but it's clearly a plastid record, coerce to Plant
    if ($kingdom_label eq 'Unknown' && $type_label eq 'Plastid') {
        $kingdom_label = 'Plant';
    }

    my $category = $kingdom_label . '-' . $type_label;

    # Track potential conflicts (same accession ? different category)
    if (exists $acc_category{$acc} && $acc_category{$acc} ne $category) {
        $conflict_count++;
        # We still skip this record as a duplicate
        $dup_count++;
        return;
    }

    $seen_acc{$acc}      = 1;
    $acc_category{$acc}  = $category;
    $cat_counts{$category}++;

    # Open output filehandle on first use
    if (!exists $fh_for_cat{$category}) {
        my $out_file = $out_prefix . '.' . $category . '.fasta';
        open my $outfh, '>>', $out_file
            or die "Cannot open output FASTA [$out_file]: $!\n";
        $fh_for_cat{$category} = $outfh;
    }

    my $fh = $fh_for_cat{$category};
    print $fh $header, "\n";
    print $fh $seq if defined $seq && $seq ne '';
}

while (my $line = <$ffh>) {
    chomp $line;
    if ($line =~ /^>/) {
        # flush previous record
        flush_record($current_header, $current_seq) if defined $current_header;
        $current_header = $line;
        $current_seq    = '';
    } else {
        # sequence line
        $current_seq .= $line . "\n";
    }
}

# flush last record
flush_record($current_header, $current_seq) if defined $current_header;

close $ffh;

# Close all output filehandles
for my $cat (keys %fh_for_cat) {
    close $fh_for_cat{$cat};
}

# Summary
print "Done splitting by kingdom + genome/marker type.\n";
print "Duplicates skipped by accession: $dup_count\n";
print "Classification conflicts (same accession, different class): $conflict_count\n";
print "Category counts:\n";
for my $cat (sort keys %cat_counts) {
    printf "  %s : %d\n", $cat, $cat_counts{$cat};
}
