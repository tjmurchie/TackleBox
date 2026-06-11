#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);

# Split a combined NCBI FASTA (from NCBI-NT_Downloader.pl) into separate files
# by:
#   - Kingdom  (Animal / Plant / Fungi / Bacteria / Archaea / Protist / Protozoa / Unknown)
#   - Type     (Mito / Plastid / NucMark / Other)
#   - Phylum   (optional; adds the phylum name to the filename when --split-phylum is used)
#
# Output filenames (default):
#   OUTPREFIX.Animal-Mito.fasta
#   OUTPREFIX.Plant-Plastid.fasta
#   OUTPREFIX.Unknown-Other.fasta
#
# Output filenames with --split-phylum:
#   OUTPREFIX.Animal.Arthropoda-Mito.fasta
#   OUTPREFIX.Plant.Tracheophyta-Plastid.fasta
#   OUTPREFIX.Unknown.Unknown-Other.fasta     (when phylum is not in the map)
#
# Usage:
#   perl split_fasta_by_kingdom_organelle_simple.pl [--split-phylum] \
#       INPUT.fasta OUTPREFIX species_kingdom.tsv [regions_config.tsv]
#
# species_kingdom.tsv format (from gbif_prep_from_csv.py):
#   species<TAB>kingdom[<TAB>phylum]

# ── CLI ──────────────────────────────────────────────────────────────────

my $split_phylum = 0;

GetOptions(
    "split-phylum" => \$split_phylum,
) or die "Usage: $0 [--split-phylum] INPUT.fasta OUTPREFIX species_kingdom.tsv [regions_config.tsv]\n";

my ($fasta_in, $out_prefix, $sp_kd_tsv, $regions_conf) = @ARGV;
if (!defined $fasta_in || !defined $out_prefix || !defined $sp_kd_tsv) {
    die "Usage: $0 [--split-phylum] INPUT.fasta OUTPREFIX species_kingdom.tsv [regions_config.tsv]\n";
}

# ── Load species -> kingdom [+phylum] map ─────────────────────────────────

my %species_to_kd;    # lc(species) -> kingdom
my %genus_to_kd;      # lc(genus)   -> kingdom  (fallback)
my %species_to_ph;    # lc(species) -> phylum   (for --split-phylum)
my %genus_to_ph;      # lc(genus)   -> phylum   (fallback)
my $rows = 0;

open my $tfh, '<', $sp_kd_tsv
    or die "Cannot open species->kingdom TSV [$sp_kd_tsv]: $!\n";

while (my $line = <$tfh>) {
    chomp $line;
    next unless $line =~ /\S/;
    my @cols = split /\t/, $line, 3;
    my ($sp, $kd, $ph) = @cols;
    next unless defined $sp && defined $kd;
    s/^\s+|\s+$//g for ($sp, $kd);
    $ph //= '';
    $ph =~ s/^\s+|\s+$//g;
    next unless $sp && $kd;

    $rows++;
    my $sp_lc = lc $sp;
    $species_to_kd{$sp_lc} = $kd;
    $species_to_ph{$sp_lc} = $ph if $ph;

    my ($genus) = split /\s+/, $sp;
    if ($genus) {
        my $g_lc = lc $genus;
        $genus_to_kd{$g_lc} //= $kd;
        $genus_to_ph{$g_lc} //= $ph if $ph;
    }
}
close $tfh;

warn sprintf("Loaded species_kingdom map: %d rows, %d species, %d genus keys%s\n",
    $rows,
    scalar(keys %species_to_kd),
    scalar(keys %genus_to_kd),
    $split_phylum ? " (phylum splitting ON)" : "");

# ── Helpers ───────────────────────────────────────────────────────────────

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

sub normalize_phylum_label {
    my ($ph_raw) = @_;
    return 'Unknown' unless defined $ph_raw && $ph_raw =~ /\S/;
    # Capitalize first letter, lowercase rest
    $ph_raw =~ s/^\s+|\s+$//g;
    $ph_raw = ucfirst(lc $ph_raw);
    # Remove characters that would break filenames
    $ph_raw =~ s{[/\\:*?"<>|]}{}g;
    return $ph_raw || 'Unknown';
}

sub classify_type {
    my ($header_lc) = @_;
    return 'Mito'    if $header_lc =~ /mitochondr/ || $header_lc =~ /\[location=mitochondrion\]/;
    return 'Plastid' if $header_lc =~ /chloroplast|plastid|plastome/
                     || $header_lc =~ /\[location=(chloroplast|plastid)\]/;
    return 'NucMark' if $header_lc =~ /\b18s\b/
                     || $header_lc =~ /\b28s\b/
                     || $header_lc =~ /small subunit ribosomal|ssu rna/
                     || $header_lc =~ /large subunit ribosomal|lsu rna/
                     || $header_lc =~ /\bits1\b|\bits2\b|internal transcribed spacer/
                     || $header_lc =~ /5\.8s/
                     || $header_lc =~ /histone h3|\bh3\b/;
    return 'Other';
}

sub parse_accession {
    my ($header) = @_;
    $header =~ s/^>//;
    my @f = split /\s+/, $header;
    return $f[0] if @f;
    return undef;
}

sub guess_species_name {
    my ($header) = @_;
    $header =~ s/^>//;
    my @f = split /\s+/, $header;
    return undef if @f < 3;
    # NCBI: ">ACC.1 Genus species description..."
    my $genus = $f[1];
    my $sp_ep = $f[2];
    return undef unless $genus =~ /^[A-Z]/;
    return "$genus $sp_ep";
}

# ── Main splitting ────────────────────────────────────────────────────────

open my $ffh, '<', $fasta_in
    or die "Cannot open FASTA file [$fasta_in]: $!\n";

my %seen_acc;
my %fh_for_cat;
my %cat_counts;
my $dup_count      = 0;
my $current_header = undef;
my $current_seq    = '';

sub flush_record {
    my ($header, $seq) = @_;
    return unless defined $header && $header ne '';

    my $acc = parse_accession($header);
    return unless $acc;

    if (exists $seen_acc{$acc}) {
        $dup_count++;
        return;
    }

    my $header_lc = lc $header;

    # Determine kingdom
    my $species_guess = guess_species_name($header);
    my ($kd_raw, $ph_raw);

    if ($species_guess) {
        my $sp_lc = lc $species_guess;
        if (exists $species_to_kd{$sp_lc}) {
            $kd_raw = $species_to_kd{$sp_lc};
            $ph_raw = $species_to_ph{$sp_lc};
        } else {
            my ($genus) = split /\s+/, $species_guess;
            if ($genus) {
                my $g_lc = lc $genus;
                $kd_raw = $genus_to_kd{$g_lc} if exists $genus_to_kd{$g_lc};
                $ph_raw = $genus_to_ph{$g_lc} if exists $genus_to_ph{$g_lc};
            }
        }
    }

    my $kingdom_label = normalize_kingdom_label($kd_raw);
    my $type_label    = classify_type($header_lc);

    # Coerce plastid-without-kingdom to Plant
    if ($kingdom_label eq 'Unknown' && $type_label eq 'Plastid') {
        $kingdom_label = 'Plant';
    }

    my $category;
    if ($split_phylum) {
        my $phylum_label = normalize_phylum_label($ph_raw);
        $category = $kingdom_label . '.' . $phylum_label . '-' . $type_label;
    } else {
        $category = $kingdom_label . '-' . $type_label;
    }

    $seen_acc{$acc} = 1;
    $cat_counts{$category}++;

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
        flush_record($current_header, $current_seq) if defined $current_header;
        $current_header = $line;
        $current_seq    = '';
    } else {
        $current_seq .= $line . "\n";
    }
}
flush_record($current_header, $current_seq) if defined $current_header;

close $ffh;
close $fh_for_cat{$_} for keys %fh_for_cat;

print "Done splitting by kingdom" . ($split_phylum ? " + phylum" : "") . " + region type.\n";
print "Duplicates skipped by accession: $dup_count\n";
print "Category counts:\n";
for my $cat (sort keys %cat_counts) {
    printf "  %-40s : %d\n", $cat, $cat_counts{$cat};
}
