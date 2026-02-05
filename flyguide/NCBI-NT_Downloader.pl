#!/usr/bin/env perl

# Download mitochondrial + plastid/chloroplast + selected nuclear marker
# NCBI nucleotide records for a list of species (one per line) and write
# them in FASTA format.
#
# For each input line:
#   1) Normalize to "Genus species" (first two tokens), unless the line
#      is a taxid like "txid12345", in which case it is used as-is.
#
#   2) Run a single organelle + nuclear-marker-focused query:
#
#        (Genus species[ORGN])
#        AND 50:400000[SLEN]
#        AND (
#          mitochondrion OR mitochondrial OR chloroplast OR plastid
#          OR (18S/28S/ITS/5.8S/H3 in [Title])
#        )
#        NOT wgs[PROP]
#        NOT tsa[PROP]
#        NOT clone[Title]
#        NOT UNVERIFIED[Title]
#        NOT chromosome[Title]
#        NOT PREDICTED[Title]
#
# Output:
#   <OUTPREFIX>.fasta         - all matching sequences in FASTA
#   <OUTPREFIX>.completed.txt - tab-separated: raw_species_line <TAB> num_records
#   <OUTPREFIX>.progress.txt  - human-readable progress + ETA
#
# Usage:
#   perl NCBI-NT_Downloader.pl species.txt OUTPREFIX you@uni.edu [NCBI_API_KEY]
#
# Requires:
#   - Perl
#   - Bio::DB::EUtilities  (e.g. via conda: perl-bio-eutilities)

use strict;
use warnings;
use Bio::DB::EUtilities;
use IO::Handle;

# ------------- Command line -------------

my $species_file = shift @ARGV or die "Usage: $0 species.txt OUTPREFIX you\@email [NCBI_API_KEY]\n";
my $out_prefix   = shift @ARGV or die "Usage: $0 species.txt OUTPREFIX you\@email [NCBI_API_KEY]\n";
my $email        = shift @ARGV || 'youremail@example.com';
my $api_key      = shift @ARGV;  # optional, but recommended

# ------------- Tuning -------------

my $SLEEP_BETWEEN_CALLS = 0.4;   # seconds between requests (~3 req/s with API key)
my $BATCH_SIZE_IDS      = 800;   # max IDs to fetch per EFETCH & per ESEARCH
my $MAX_RETRIES         = 3;     # retry transient errors this many times

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
    $line =~ s/^\s+//;
    $line =~ s/\s+$//;
    next unless $line;
    push @species_raw, $line;
}
close $sfh;

# --- progress tracking ---
my $total_species  = scalar @species_raw;
my $already_done   = scalar keys %done;
my $start_time     = time;
my $processed_now  = 0;

# ------------- Open FASTA output -------------

open my $outfh, '>>', "$out_prefix.fasta"
    or die "Cannot open $out_prefix.fasta for writing: $!\n";

# ------------- Main loop -------------

for my $raw_species (@species_raw) {

    if ($done{$raw_species}) {
        warn "Skipping already completed species line: $raw_species\n";
        next;
    }

    # Normalize GBIF-style names: keep only "Genus species"
    my $norm_species = normalize_species_name($raw_species);

    warn "=== Processing [$raw_species] as organism [$norm_species] ===\n";

    # Build ORGN field term: e.g. Sanogasta backhauseni[ORGN]
    my $org_term = $norm_species . "[ORGN]";

    # Organelle + nuclear markers query
    my $term =
        "($org_term) " .
        "AND 50:400000[SLEN] " .
        "AND (" .
          # 1) Organelle sequences
          "(mitochondrion OR mitochondrial OR chloroplast OR plastid) " .
          # 2) OR nuclear markers (18S/28S/ITS/H3) in the title
          "OR (" .
            "(" .
              "\"18S ribosomal RNA\"[Title] OR 18S[Title] OR " .
              "\"small subunit ribosomal\"[Title] OR \"SSU rRNA\"[Title] OR " .
              "\"28S ribosomal RNA\"[Title] OR 28S[Title] OR " .
              "\"large subunit ribosomal\"[Title] OR \"LSU rRNA\"[Title] OR " .
              "ITS1[Title] OR ITS2[Title] OR \"internal transcribed spacer\"[Title] OR " .
              "\"5.8S ribosomal RNA\"[Title] OR \"5.8S rRNA\"[Title] OR " .
              "\"histone H3\"[Title] OR H3[Title]" .
            ") " .
          ")" .
        ") " .
        "NOT wgs[PROP] " .
        "NOT tsa[PROP] " .
        "NOT clone[Title] " .
        "NOT UNVERIFIED[Title] " .
        "NOT chromosome[Title] " .
        "NOT PREDICTED[Title]";

    warn "ESearch term: $term\n";

    my ($count, @ids) = esearch_ids(
        $term,
        $email, $api_key,
        $BATCH_SIZE_IDS,
        $SLEEP_BETWEEN_CALLS,
        $MAX_RETRIES
    );
    $count //= 0;
    warn "ESearch count for [$norm_species]: $count\n";

    if (!@ids) {
        warn "No organelle/nuclear marker records for [$raw_species] (normalized: [$norm_species])\n";
        print $logfh "$raw_species\t0\n";

        # --- progress update even for 0-record taxa ---
        $processed_now++;

        my $completed = $already_done + $processed_now;
        my $elapsed   = time - $start_time;
        my $frac      = $total_species ? $completed / $total_species : 0;
        my $remaining = $total_species - $completed;
        $remaining = 0 if $remaining < 0;

        my $eta_secs  = ($processed_now > 0 && $remaining > 0)
            ? int( $elapsed / $processed_now * $remaining )
            : 0;

        my $elapsed_str = format_time($elapsed);
        my $eta_str     = format_time($eta_secs);

        warn sprintf("Progress: %d/%d (%.1f%%). Elapsed: %s, ETA: %s\n",
                     $completed, $total_species, $frac * 100,
                     $elapsed_str, $eta_str);

        if (open my $pfh, '>', "$out_prefix.progress.txt") {
            printf $pfh "Completed: %d/%d (%.1f%%)\nElapsed: %s\nETA: %s\n",
                $completed, $total_species, $frac * 100,
                $elapsed_str, $eta_str;
            close $pfh;
        }

        sleep $SLEEP_BETWEEN_CALLS;
        next;
    }

    my $n = scalar @ids;
    warn "Fetching $n organelle/nuclear marker records for [$raw_species]\n";

    efetch_and_write(
        \@ids,
        $outfh,
        $email, $api_key,
        $BATCH_SIZE_IDS,
        $SLEEP_BETWEEN_CALLS,
        $MAX_RETRIES
    );

    print $logfh "$raw_species\t$n\n";

    # --- progress update after each species with hits ---
    $processed_now++;

    my $completed = $already_done + $processed_now;
    my $elapsed   = time - $start_time;
    my $frac      = $total_species ? $completed / $total_species : 0;
    my $remaining = $total_species - $completed;
    $remaining = 0 if $remaining < 0;

    my $eta_secs  = ($processed_now > 0 && $remaining > 0)
        ? int( $elapsed / $processed_now * $remaining )
        : 0;

    my $elapsed_str = format_time($elapsed);
    my $eta_str     = format_time($eta_secs);

    warn sprintf("Progress: %d/%d (%.1f%%). Elapsed: %s, ETA: %s\n",
                 $completed, $total_species, $frac * 100,
                 $elapsed_str, $eta_str);

    if (open my $pfh, '>', "$out_prefix.progress.txt") {
        printf $pfh "Completed: %d/%d (%.1f%%)\nElapsed: %s\nETA: %s\n",
            $completed, $total_species, $frac * 100,
            $elapsed_str, $eta_str;
        close $pfh;
    }

    sleep $SLEEP_BETWEEN_CALLS;
}

close $outfh;
close $logfh;

exit 0;

# ------------- Subroutines -------------

# Normalize a GBIF-style species string to "Genus species"
# e.g. "Sanogasta backhauseni (Simon, 1895)" -> "Sanogasta backhauseni"
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

# Return (count, @ids) for a given term
sub esearch_ids {
    my ($term, $email, $api_key, $retmax, $sleep, $max_retries) = @_;

    my @ids;
    my $count   = 0;
    my $attempt = 0;

    while ($attempt < $max_retries) {
        $attempt++;
        eval {
            my %params = (
                -eutil  => 'esearch',
                -db     => 'nucleotide',
                -term   => $term,
                -email  => $email,
                -retmax => $retmax,
            );
            $params{'-api_key'} = $api_key if defined $api_key && length $api_key;

            my $fac = Bio::DB::EUtilities->new(%params);
            $count = $fac->get_count;
            @ids   = $fac->get_ids;
        };
        if ($@) {
            warn "ESearch error (attempt $attempt/$max_retries) for term [$term]: $@\n";
            sleep($sleep * $attempt);
        } else {
            last;
        }
    }

    return ($count, @ids);
}

sub efetch_and_write {
    my ($uids_ref, $outfh, $email, $api_key,
        $chunk_size, $sleep, $max_retries) = @_;

    my @uids = @$uids_ref;

    while (@uids) {
        my @chunk = splice(@uids, 0, $chunk_size);

        my $attempt = 0;
        my $content;

        while ($attempt < $max_retries) {
            $attempt++;
            eval {
                my %params = (
                    -eutil   => 'efetch',
                    -db      => 'nucleotide',
                    -id      => \@chunk,
                    -rettype => 'fasta',
                    -retmode => 'text',
                    -email   => $email,
                );
                $params{'-api_key'} = $api_key if defined $api_key && length $api_key;

                my $fac = Bio::DB::EUtilities->new(%params);
                $content = $fac->get_Response->content;
            };
            if ($@) {
                warn "EFetch error (attempt $attempt/$max_retries) for "
                     . scalar(@chunk) . " IDs: $@\n";
                sleep($sleep * $attempt);
            } else {
                last;
            }
        }

        if (defined $content) {
            print $outfh $content;
        } else {
            warn "Giving up on a chunk of " . scalar(@chunk) .
                 " IDs after repeated efetch failures\n";
        }

        sleep $sleep;
    }
}

sub format_time {
    my ($secs) = @_;
    $secs = int($secs);
    my $h = int($secs / 3600);
    $secs -= $h * 3600;
    my $m = int($secs / 60);
    $secs -= $m * 60;
    return sprintf("%02dh:%02dm:%02ds", $h, $m, $secs);
}
