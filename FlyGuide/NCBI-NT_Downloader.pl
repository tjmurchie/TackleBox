#!/usr/bin/env perl

# TackleBox: FlyGuide - NCBI-NT_Downloader
#
# Download mitochondrial + plastid/chloroplast + selected marker
# NCBI nucleotide records for a list of species (one per line) and write
# them in FASTA format.
#
# For each input line:
#   1) Normalize to "Genus species" (first two tokens), unless the line
#      is a taxid like "txid12345", in which case it is used as-is.
#
#   2) Build a query of the form:
#
#        (Genus species[ORGN])
#        AND 50:400000[SLEN]
#        AND (
#          (mitochondrion OR mitochondrial OR chloroplast OR plastid)
#          OR <marker Title clauses>
#        )
#        NOT wgs[PROP]
#        NOT tsa[PROP]
#        NOT clone[Title]
#        NOT UNVERIFIED[Title]
#        NOT chromosome[Title]
#        NOT PREDICTED[Title]
#
#      - If regions_config.tsv is present and usable, <marker Title clauses>
#        are built from enabled rows with a non-empty ncbi_title_clause.
#      - Otherwise, we fall back to a hard-coded default:
#        18S / 28S / ITS / 5.8S / histone H3 Title clauses.
#
#   3) ESearch to get matching IDs (up to --max-per-taxon).
#   4) EFETCH in batches and append FASTA records to OUTPREFIX.fasta.
#
# CLI:
#   NCBI-NT_Downloader.pl [--max-per-taxon N] species.txt OUTPREFIX you@email [NCBI_API_KEY]
#
#   --max-per-taxon N   Maximum number of records to fetch per taxon (default: 1000)

use strict;
use warnings;
use Bio::DB::EUtilities;
use IO::Handle;
use FindBin;
use Getopt::Long qw(GetOptions);

# Optional configuration for Title clauses (markers) via regions_config.tsv
my $REGIONS_CONFIG = "$FindBin::Bin/regions_config.tsv";
my @MARKER_CLAUSES;         # NCBI [Title] clauses from regions_config.tsv
my $USE_CONFIG_MARKERS = 0; # set to 1 if config successfully loaded

# ------------- CLI options -------------

my $max_per_taxon = 1000;   # default limit per taxon

GetOptions(
    "max-per-taxon=i" => \$max_per_taxon,
) or die "Usage: $0 [--max-per-taxon N] species.txt OUTPREFIX you\@email [NCBI_API_KEY]\n";

# ------------- Positional args -------------

my $species_file = shift @ARGV
  or die "Usage: $0 [--max-per-taxon N] species.txt OUTPREFIX you\@email [NCBI_API_KEY]\n";
my $out_prefix   = shift @ARGV
  or die "Usage: $0 [--max-per-taxon N] species.txt OUTPREFIX you\@email [NCBI_API_KEY]\n";
my $email        = shift @ARGV || 'youremail@example.com';
my $api_key      = shift @ARGV;  # optional but recommended

# ------------- Tuning -------------

my $SLEEP_BETWEEN_CALLS = 0.4;   # seconds between requests (~3 req/s with API key)
my $MAX_RETRIES         = 3;

# ------------- Load marker clauses (config + fallback) -------------

load_marker_clauses_from_regions($REGIONS_CONFIG);

# If config didn't give us any clauses, use a hard-coded default nuclear
# marker block (18S/28S/ITS/5.8S/H3).
my $DEFAULT_MARKER_BLOCK = join(" OR ", (
    "\"18S ribosomal RNA\"[Title]",
    "18S[Title]",
    "\"small subunit ribosomal\"[Title]",
    "\"SSU rRNA\"[Title]",
    "\"28S ribosomal RNA\"[Title]",
    "28S[Title]",
    "\"large subunit ribosomal\"[Title]",
    "\"LSU rRNA\"[Title]",
    "ITS1[Title]",
    "ITS2[Title]",
    "\"internal transcribed spacer\"[Title]",
    "\"5.8S ribosomal RNA\"[Title]",
    "\"5.8S rRNA\"[Title]",
    "\"histone H3\"[Title]",
    "H3[Title]"
));

if (!$USE_CONFIG_MARKERS) {
    warn "[NCBI-NT_Downloader] No usable Title clauses from regions_config.tsv; using built-in nuclear marker block (18S/28S/ITS/5.8S/H3).\n";
}

# ------------- Restart bookkeeping -------------

my %done;
my $done_file = $out_prefix . ".completed.txt";

if (-e $done_file) {
    open my $dfh, '<', $done_file or die "Cannot open $done_file for reading: $!\n";
    while (my $line = <$dfh>) {
        chomp $line;
        next unless $line =~ /\S/;
        my ($sp_raw) = split /\t/, $line;
        $done{$sp_raw} = 1;
    }
    close $dfh;
}

open my $logfh, '>>', $done_file or die "Cannot open $done_file for appending: $!\n";
$logfh->autoflush(1);

# ------------- Read species list -------------

open my $sfh, '<', $species_file or die "Cannot open $species_file: $!\n";
my @species_raw;
while (my $line = <$sfh>) {
    chomp $line;
    next unless $line =~ /\S/;
    push @species_raw, $line;
}
close $sfh;

my $total_species = scalar @species_raw;
warn "Loaded $total_species species lines from $species_file\n";
warn "Max records per taxon (ESearch/EFetch): $max_per_taxon\n";

my $already_done = scalar keys %done;
if ($already_done) {
    warn "Found $already_done completed lines in $done_file; will skip those.\n";
}

my $start_time    = time;
my $processed_now = 0;

# ------------- Open FASTA output -------------

open my $outfh, '>>', "$out_prefix.fasta"
    or die "Cannot open $out_prefix.fasta for writing: $!\n";

# ------------- Main loop -------------

for my $raw_species (@species_raw) {

    if ($done{$raw_species}) {
        warn "Skipping already completed species line: $raw_species\n";
        next;
    }

    my $norm_species = normalize_species_name($raw_species);
    my $org_term;

    if ($norm_species =~ /^txid\d+$/i) {
        # Direct taxid line like "txid2172571"
        $org_term = $norm_species . "[ORGN]";
    } else {
        $org_term = $norm_species . "[ORGN]";
    }

    # Organelle block (mitochondrial + plastid)
    my $organelle_block = "(mitochondrion OR mitochondrial OR chloroplast OR plastid)";

    my $term;

    if ($USE_CONFIG_MARKERS && @MARKER_CLAUSES) {
        # Config-driven behaviour: organelle block plus any enabled Title clauses
        my $marker_block = "(" . join(" OR ", @MARKER_CLAUSES) . ")";

        $term =
            "($org_term) " .
            "AND 50:400000[SLEN] " .
            "AND (" .
              $organelle_block . " OR " . $marker_block .
            ") " .
            "NOT wgs[PROP] " .
            "NOT tsa[PROP] " .
            "NOT clone[Title] " .
            "NOT UNVERIFIED[Title] " .
            "NOT chromosome[Title] " .
            "NOT PREDICTED[Title]";
    } else {
        # Full default behaviour: organelle + hard-coded nuclear marker block
        $term =
            "($org_term) " .
            "AND 50:400000[SLEN] " .
            "AND (" .
              $organelle_block . " OR (" . $DEFAULT_MARKER_BLOCK . ")" .
            ") " .
            "NOT wgs[PROP] " .
            "NOT tsa[PROP] " .
            "NOT clone[Title] " .
            "NOT UNVERIFIED[Title] " .
            "NOT chromosome[Title] " .
            "NOT PREDICTED[Title]";
    }

    warn "ESearch term: $term\n";

    my ($count, @ids) = esearch_ids(
        $term,
        $email, $api_key,
        $max_per_taxon,          # retmax = max records per taxon
        $SLEEP_BETWEEN_CALLS,
        $MAX_RETRIES
    );
    $count //= 0;
    warn "ESearch count for [$norm_species]: $count\n";
    warn "  (retrieving up to $max_per_taxon IDs for this taxon)\n";

    if (!@ids) {
        warn "No organelle/marker records for [$raw_species] (normalized: [$norm_species])\n";
        print $logfh "$raw_species\t0\n";

        $processed_now++;
        update_progress($total_species, $already_done, $processed_now, $start_time, $out_prefix);

        sleep $SLEEP_BETWEEN_CALLS;
        next;
    }

    # Fetch in batches (chunk size = max_per_taxon, so usually 1 batch)
    my $chunk_size = $max_per_taxon;
    my $batches = int((scalar(@ids) + $chunk_size - 1) / $chunk_size);
    warn "Fetching in $batches batch(es)\n";

    for (my $i = 0; $i < @ids; $i += $chunk_size) {

        my @chunk = @ids[$i .. ($i + $chunk_size - 1 < $#ids ? $i + $chunk_size - 1 : $#ids)];
        my $batch_index = int($i / $chunk_size) + 1;

        warn "  EFETCH batch $batch_index / $batches (", scalar(@chunk), " IDs)\n";

        my $fa = efetch_batch(
            \@chunk,
            $email, $api_key,
            $SLEEP_BETWEEN_CALLS,
            $MAX_RETRIES
        );

        if (defined $fa && $fa ne '') {
            print $outfh $fa;
        } else {
            warn "  EFETCH batch $batch_index returned no data\n";
        }
    }

    # Mark this species line as completed
    print $logfh "$raw_species\t$count\n";

    $processed_now++;
    update_progress($total_species, $already_done, $processed_now, $start_time, $out_prefix);

    sleep $SLEEP_BETWEEN_CALLS;
}

close $outfh;
close $logfh;

exit 0;

# ------------- Subroutines -------------

# Load Title clauses (for NCBI [Title] queries) from regions_config.tsv.
# Any row where enabled_default != 0 and ncbi_title_clause is non-empty
# contributes a clause (regardless of class).
sub load_marker_clauses_from_regions {
    my ($path) = @_;
    return unless defined $path && -f $path;

    open my $rfh, '<', $path or do {
        warn "[NCBI-NT_Downloader] Could not open regions config [$path]: $!; using built-in nuclear marker block.\n";
        return;
    };

    my $line_no = 0;
    my @clauses;
    while (my $line = <$rfh>) {
        chomp $line;
        next if $line =~ /^\s*#/;      # comment
        next unless $line =~ /\S/;     # blank

        my @cols = split /\t/, $line, 5;
        # Header row?
        if ($line_no == 0 && $cols[0] =~ /region_id/i) {
            $line_no++;
            next;
        }
        $line_no++;

        my ($id, $class, $enabled, $regex, $title_clause) = @cols;
        $enabled      //= 1;
        $title_clause //= '';

        next if $enabled eq '0';
        next unless $title_clause =~ /\S/;

        $title_clause =~ s/^\s+//;
        $title_clause =~ s/\s+$//;

        push @clauses, "(" . $title_clause . ")";
    }
    close $rfh;

    if (@clauses) {
        @MARKER_CLAUSES     = @clauses;
        $USE_CONFIG_MARKERS = 1;
        warn "[NCBI-NT_Downloader] Loaded " . scalar(@MARKER_CLAUSES) . " Title clauses from $path\n";
    } else {
        $USE_CONFIG_MARKERS = 0;
    }
}

# Normalize a species string to "Genus species", unless it's a txid.
sub normalize_species_name {
    my ($raw) = @_;

    $raw =~ s/^\s+//;
    $raw =~ s/\s+$//;

    # Allow direct taxid lines like "txid2172571"
    if ($raw =~ /^txid\d+$/i) {
        return $raw;
    }

    my @parts = split(/\s+/, $raw);
    if (@parts >= 2) {
        return $parts[0] . " " . $parts[1];
    }
    return $raw;
}

# Return (count, @ids) for a given term.
sub esearch_ids {
    my ($term, $email, $api_key, $retmax, $sleep, $max_retries) = @_;

    my @ids;
    my $count   = 0;
    my $attempt = 0;

    while ($attempt < $max_retries) {
        $attempt++;

        my $eutil = Bio::DB::EUtilities->new(
            -eutil   => 'esearch',
            -db      => 'nuccore',
            -term    => $term,
            -email   => $email,
            -api_key => $api_key,
            -retmax  => $retmax,
        );

        eval {
            $count = $eutil->get_count;
            @ids   = $eutil->get_ids;
        };
        if ($@) {
            warn "ESearch attempt $attempt failed: $@\n";
            if ($attempt < $max_retries) {
                warn "Sleeping $sleep seconds before retry...\n";
                sleep $sleep;
                next;
            } else {
                warn "Max retries reached for ESearch; giving up.\n";
                last;
            }
        }

        last;
    }

    return ($count, @ids);
}

# EFETCH a batch of IDs and return FASTA text.
sub efetch_batch {
    my ($id_list, $email, $api_key, $sleep, $max_retries) = @_;

    my $attempt = 0;
    my $fasta  = '';

    while ($attempt < $max_retries) {
        $attempt++;

        my $eutil = Bio::DB::EUtilities->new(
            -eutil   => 'efetch',
            -db      => 'nuccore',
            -rettype => 'fasta',
            -retmode => 'text',
            -id      => $id_list,
            -email   => $email,
            -api_key => $api_key,
        );

        eval {
            $fasta = $eutil->get_Response->content;
        };
        if ($@) {
            warn "EFetch attempt $attempt failed: $@\n";
            if ($attempt < $max_retries) {
                warn "Sleeping $sleep seconds before retry...\n";
                sleep $sleep;
                next;
            } else {
                warn "Max retries reached for EFetch; giving up on this batch of IDs.\n";
                last;
            }
        }

        last;
    }

    return $fasta;
}

# Progress reporting & ETA
sub update_progress {
    my ($total_species, $already_done, $processed_now, $start_time, $prefix) = @_;

    my $completed = $already_done + $processed_now;
    my $elapsed   = time - $start_time;
    my $frac      = $total_species ? $completed / $total_species : 0;

    my $elapsed_str = format_time($elapsed);
    my $eta_seconds = $frac > 0 ? $elapsed / $frac - $elapsed : 0;
    my $eta_str     = format_time($eta_seconds);

    warn sprintf("Progress: %d/%d (%.1f%%); elapsed: %s; ETA: %s\n",
                 $completed, $total_species, $frac * 100,
                 $elapsed_str, $eta_str);

    if (open my $pfh, '>', "$prefix.progress.txt") {
        printf $pfh "Completed: %d/%d (%.1f%%)\nElapsed: %s\nETA: %s\n",
            $completed, $total_species, $frac * 100,
            $elapsed_str, $eta_str;
        close $pfh;
    }
}

sub format_time {
    my ($seconds) = @_;
    $seconds = int($seconds);
    my $h = int($seconds / 3600);
    my $m = int(($seconds % 3600) / 60);
    my $s = $seconds % 60;
    return sprintf("%02dh:%02dm:%02ds", $h, $m, $s);
}
