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
#        (Genus species[ORGN])
#        AND MIN:MAX[SLEN]
#        AND (organelle_block OR marker_block)
#        NOT wgs[PROP] NOT tsa[PROP] NOT clone[Title]
#        NOT UNVERIFIED[Title] NOT chromosome[Title] NOT PREDICTED[Title]
#
#      Marker block is built from regions_config.tsv (ncbi_title_clause column)
#      or falls back to built-in 18S/28S/ITS/5.8S/H3 clauses.
#
#   3) ESearch to get matching IDs (up to --max-per-taxon).
#   4) EFETCH in batches and append FASTA records to OUTPREFIX.fasta.
#
# CLI:
#   NCBI-NT_Downloader.pl [OPTIONS] species.txt OUTPREFIX you@email [NCBI_API_KEY]
#
# Options:
#   --max-per-taxon N      Max NCBI records to fetch per taxon (default: 1000)
#   --min-slen N           Min sequence length filter in bp (default: 50)
#   --max-slen N           Max sequence length filter in bp (default: 400000)
#   --no-organelle         Skip organelle records; fetch nuclear markers only
#   --enable-regions LIST  Comma-separated region_id(s) to force-enable
#                          (overrides enabled_default=0 in regions_config.tsv)
#   --disable-regions LIST Comma-separated region_id(s) to force-disable
#   --no-tui               Disable the interactive TUI; use plain line output
#   --name-mode MODE        Taxon normalization mode for [ORGN] query:
#                             species   = legacy first two tokens (default)
#                             trinomial = first three tokens if present
#                             as-is     = full input line after whitespace cleanup
#
# Query customization via regions_config.tsv:
#   Each row in regions_config.tsv contributes an NCBI [Title] search clause.
#   Rows with enabled_default=0 are skipped unless --enable-regions is used.
#   Example: --enable-regions COI        (adds COI/cox1 to every query)
#            --disable-regions NUC_H3    (removes histone H3 from queries)
#            --min-slen 100 --max-slen 5000  (tighter length filter)
#            --no-organelle              (nuclear markers only)

use strict;
use warnings;
use Bio::DB::EUtilities;
use IO::Handle;
use FindBin;
use Getopt::Long qw(GetOptions);

STDERR->autoflush(1);
STDOUT->autoflush(1);
binmode(STDERR, ':encoding(UTF-8)');

# ── Terminal detection ────────────────────────────────────────────────────
my $IS_TTY = -t STDERR ? 1 : 0;

my $TERM_WIDTH = 80;
if ($IS_TTY) {
    my $w = `tput cols 2>/dev/null`;
    $w //= '';
    chomp($w);
    $TERM_WIDTH = ($w =~ /^\d+$/ && $w > 30) ? int($w) : 80;
}

# Dashboard state (file-scope so subroutines can access)
my $RECENT_MAX    = 5;
my @recent_taxa   = ();   # array of { name=>, count=> } hashrefs
my $current_taxon  = '';
my $current_status = '';
my $dash_drawn     = 0;
my $burnin         = 10;  # recalculated after loading species list

# ── Signal handler: restore terminal on Ctrl+C ───────────────────────────
$SIG{INT} = sub {
    if ($IS_TTY) {
        print STDERR "\033[?25h\n";   # restore cursor visibility
    }
    exit 130;
};

# ── Regions configuration ─────────────────────────────────────────────────
my $REGIONS_CONFIG   = "$FindBin::Bin/regions_config.tsv";
my @MARKER_CLAUSES;
my $USE_CONFIG_MARKERS = 0;

# ── CLI options ───────────────────────────────────────────────────────────
my $max_per_taxon     = 1000;
my $min_slen          = 50;
my $max_slen          = 400000;
my $no_organelle      = 0;
my $no_tui            = 0;
my $enable_regions_s  = '';
my $disable_regions_s = '';
my $name_mode         = 'species';

GetOptions(
    "max-per-taxon=i"    => \$max_per_taxon,
    "min-slen=i"         => \$min_slen,
    "max-slen=i"         => \$max_slen,
    "no-organelle"       => \$no_organelle,
    "enable-regions=s"   => \$enable_regions_s,
    "disable-regions=s"  => \$disable_regions_s,
    "no-tui"             => \$no_tui,
    "name-mode=s"        => \$name_mode,
) or die "Usage: $0 [options] species.txt OUTPREFIX you\@email [NCBI_API_KEY]\n";

$name_mode = lc($name_mode // 'species');
die "ERROR: --name-mode must be one of: species, trinomial, as-is\n"
    unless $name_mode =~ /^(species|trinomial|as-is)$/;

$IS_TTY = 0 if $no_tui;

my %force_enable  = map { $_ => 1 } grep { /\S/ } split /,/, $enable_regions_s;
my %force_disable = map { $_ => 1 } grep { /\S/ } split /,/, $disable_regions_s;

# ── Positional args ───────────────────────────────────────────────────────
my $species_file = shift @ARGV
    or die "Usage: $0 [options] species.txt OUTPREFIX you\@email [NCBI_API_KEY]\n";
my $out_prefix = shift @ARGV
    or die "Usage: $0 [options] species.txt OUTPREFIX you\@email [NCBI_API_KEY]\n";
my $email   = shift @ARGV || 'youremail@example.com';
my $api_key = shift @ARGV;

# ── Tuning ────────────────────────────────────────────────────────────────
my $SLEEP_BETWEEN_CALLS = 0.4;
my $MAX_RETRIES         = 3;

# ── Load marker clauses from regions_config.tsv ───────────────────────────
load_marker_clauses_from_regions($REGIONS_CONFIG);

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

# ── Restart bookkeeping ───────────────────────────────────────────────────
my %done;
my $done_file = $out_prefix . ".completed.txt";

if (-e $done_file) {
    open my $dfh, '<', $done_file or die "Cannot open $done_file: $!\n";
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

# ── Read species list ─────────────────────────────────────────────────────
open my $sfh, '<', $species_file or die "Cannot open $species_file: $!\n";
my @species_raw;
while (my $line = <$sfh>) {
    chomp $line;
    next unless $line =~ /\S/;
    push @species_raw, $line;
}
close $sfh;

my $total_species = scalar @species_raw;
my $already_done  = scalar keys %done;

# ETA burn-in: show "Calculating..." until this many NEW taxa are processed.
# Uses ~2% of total, clamped to 5–20.
$burnin = int($total_species * 0.02);
$burnin = 5  if $burnin < 5;
$burnin = 20 if $burnin > 20;

# ── Print static header (once) ────────────────────────────────────────────
print_static_header();

# ── Open FASTA output ─────────────────────────────────────────────────────
open my $outfh, '>>', "$out_prefix.fasta"
    or die "Cannot open $out_prefix.fasta for writing: $!\n";

my $start_time    = time;
my $processed_now = 0;

# ── Main download loop ────────────────────────────────────────────────────
for my $raw_species (@species_raw) {

    if ($done{$raw_species}) {
        next;
    }

    my $norm_species = normalize_species_name($raw_species);
    my $org_term     = "${norm_species}[ORGN]";

    my $organelle_block = "(mitochondrion OR mitochondrial OR chloroplast OR plastid)";

    my $term;
    my $slen_filter = "${min_slen}:${max_slen}[SLEN]";
    my $not_block   = "NOT wgs[PROP] NOT tsa[PROP] NOT clone[Title] NOT UNVERIFIED[Title] NOT chromosome[Title] NOT PREDICTED[Title]";

    if ($no_organelle) {
        my $marker_block = $USE_CONFIG_MARKERS && @MARKER_CLAUSES
            ? "(" . join(" OR ", @MARKER_CLAUSES) . ")"
            : "(" . $DEFAULT_MARKER_BLOCK . ")";
        $term = "($org_term) AND $slen_filter AND $marker_block $not_block";
    } elsif ($USE_CONFIG_MARKERS && @MARKER_CLAUSES) {
        my $marker_block = "(" . join(" OR ", @MARKER_CLAUSES) . ")";
        $term = "($org_term) AND $slen_filter AND ($organelle_block OR $marker_block) $not_block";
    } else {
        $term = "($org_term) AND $slen_filter AND ($organelle_block OR ($DEFAULT_MARKER_BLOCK)) $not_block";
    }

    # Update TUI: taxon now searching
    $current_taxon  = $norm_species;
    $current_status = 'searching...';
    draw_dashboard($total_species, $already_done, $processed_now, $start_time);

    my ($count, @ids);
    {
        local $SIG{__WARN__} = $IS_TTY ? sub {} : ($SIG{__WARN__} // sub { warn @_ });
        ($count, @ids) = esearch_ids(
            $term, $email, $api_key,
            $max_per_taxon, $SLEEP_BETWEEN_CALLS, $MAX_RETRIES
        );
    }
    $count //= 0;

    if (!@ids) {
        push_recent($norm_species, 0);
        print $logfh "$raw_species\t0\n";
        $processed_now++;
        draw_dashboard($total_species, $already_done, $processed_now, $start_time);
        sleep $SLEEP_BETWEEN_CALLS;
        next;
    }

    # Fetch in batches
    my $chunk_size = $max_per_taxon;
    my $batches    = int((scalar(@ids) + $chunk_size - 1) / $chunk_size);

    for (my $i = 0; $i < @ids; $i += $chunk_size) {
        my $end   = ($i + $chunk_size - 1 < $#ids) ? $i + $chunk_size - 1 : $#ids;
        my @chunk = @ids[$i .. $end];
        my $batch_index = int($i / $chunk_size) + 1;

        if ($batches > 1) {
            $current_status = "fetching $batch_index/$batches...";
            draw_dashboard($total_species, $already_done, $processed_now, $start_time);
        }

        my $fa;
        {
            local $SIG{__WARN__} = $IS_TTY ? sub {} : ($SIG{__WARN__} // sub { warn @_ });
            $fa = efetch_batch(
                \@chunk, $email, $api_key,
                $SLEEP_BETWEEN_CALLS, $MAX_RETRIES
            );
        }
        print $outfh $fa if defined $fa && $fa ne '';
    }

    print $logfh "$raw_species\t$count\n";
    push_recent($norm_species, $count);
    $processed_now++;
    draw_dashboard($total_species, $already_done, $processed_now, $start_time);

    sleep $SLEEP_BETWEEN_CALLS;
}

close $outfh;
close $logfh;

# Final display
$current_taxon  = '';
$current_status = '';
draw_dashboard($total_species, $already_done, $processed_now, $start_time);

if ($IS_TTY) {
    print STDERR "\033[?25h\n";  # restore cursor
}

exit 0;

# ═════════════════════════════════════════════════════════════════════════════
# DISPLAY SUBROUTINES
# ═════════════════════════════════════════════════════════════════════════════

sub print_static_header {
    my $w   = $TERM_WIDTH;
    my $sep = "\x{2500}" x $w;   # ─────

    if ($IS_TTY) {
        print STDERR "\033[?25l";  # hide cursor while running
    }

    print STDERR "$sep\n";
    print STDERR "  TackleBox: FlyGuide \x{2014} NCBI Nucleotide Downloader\n";
    print STDERR "$sep\n";

    # Show query template (what SPECIES will be replaced with each taxon name)
    my $qt = build_query_template();
    my $label  = "  Query template: ";
    my $indent = " " x length($label);
    my $avail  = $w - length($label);
    if (length($qt) <= $avail) {
        print STDERR "$label$qt\n";
    } else {
        print STDERR $label . substr($qt, 0, $avail) . "\n";
        my $rest = substr($qt, $avail);
        while (length($rest) > 0) {
            my $chunk = substr($rest, 0, $w - length($indent));
            print STDERR "$indent$chunk\n";
            $rest = substr($rest, length($chunk));
        }
    }

    # Active region IDs
    my $region_info = $USE_CONFIG_MARKERS
        ? "  Markers: " . scalar(@MARKER_CLAUSES) . " clauses from regions_config.tsv"
        : "  Markers: built-in (18S/28S/ITS/5.8S/H3)";
    if (%force_enable)  { $region_info .= "  [+enabled: "  . join(",", sort keys %force_enable)  . "]"; }
    if (%force_disable) { $region_info .= "  [-disabled: " . join(",", sort keys %force_disable) . "]"; }
    print STDERR "$region_info\n";

    if ($already_done > 0) {
        print STDERR "  Resuming: $already_done / $total_species taxa already completed (will skip).\n";
    }
    print STDERR "$sep\n\n";
}

# draw_dashboard: called after every taxon. In TUI mode, moves cursor up and
# redraws a fixed-height dashboard in-place.
sub draw_dashboard {
    my ($total, $already, $processed_now, $start_time) = @_;

    my $completed = $already + $processed_now;
    my $remaining = $total - $completed;
    my $elapsed   = time - $start_time;
    my $frac      = $total > 0 ? $completed / $total : 0;
    my $w         = $TERM_WIDTH;
    my $sep       = "\x{2500}" x $w;

    my $elapsed_str = format_time($elapsed);
    my $eta_seconds = ($frac > 0) ? ($elapsed / $frac) - $elapsed : 0;
    my $eta_str = ($processed_now < $burnin)
        ? "Calculating..."
        : format_time($eta_seconds);

    # Progress bar (fills available width)
    my $pct_label = sprintf(" %5.1f%%", $frac * 100);
    my $bar_width = $w - length($pct_label) - 4;  # "  [" (3) + "]" (1) + pct_label (7) = 11 fixed
    $bar_width = 20 if $bar_width < 20;
    my $filled = int($frac * $bar_width);
    my $bar = "  ["
            . ("\x{2588}" x $filled)
            . ("\x{2591}" x ($bar_width - $filled))
            . "]" . $pct_label;

    # Build the dashboard lines
    my @lines;
    push @lines, $sep;
    push @lines, "  Recent downloads:";

    # Pad recent list to RECENT_MAX rows so dashboard height is constant
    my @display = @recent_taxa;
    unshift @display, undef while scalar(@display) < $RECENT_MAX;

    for my $entry (@display) {
        if (!defined $entry) {
            push @lines, "";
            next;
        }
        my $name  = $entry->{name};
        my $count = $entry->{count};
        my $mark  = "  \x{2713} ";  # ✓
        my $right = sprintf("[%4d records]", $count);
        my $name_w = $w - length($mark) - length($right) - 4;
        $name_w = 15 if $name_w < 15;
        $name = substr($name, 0, $name_w) if length($name) > $name_w;
        push @lines, sprintf("%s%-${name_w}s  %s", $mark, $name, $right);
    }

    # Current taxon line
    if ($current_taxon ne '') {
        my $mark  = "  \x{25B6} ";  # ▶
        my $right = sprintf("[%-14s]", $current_status);
        my $name_w = $w - length($mark) - length($right) - 4;
        $name_w = 15 if $name_w < 15;
        my $name = substr($current_taxon, 0, $name_w);
        push @lines, sprintf("%s%-${name_w}s  %s", $mark, $name, $right);
    } else {
        push @lines, "";
    }

    push @lines, $sep;
    push @lines, sprintf("  Progress: %d / %d    Remaining: %d", $completed, $total, $remaining);
    push @lines, $bar;
    push @lines, sprintf("  Elapsed: %-16s ETA: %s", $elapsed_str, $eta_str);
    push @lines, $sep;

    if ($IS_TTY) {
        if ($dash_drawn) {
            # Move cursor up to start of dashboard
            print STDERR "\033[" . scalar(@lines) . "A";
        }
        for my $line (@lines) {
            print STDERR "\033[2K$line\n";
        }
        $dash_drawn = 1;

        # Also update progress file silently
        _write_progress_file($completed, $total, $frac, $elapsed_str, $eta_str);
    } else {
        # Plain mode: one progress line per taxon
        warn sprintf("Progress: %d/%d (%.1f%%); elapsed: %s; ETA: %s\n",
            $completed, $total, $frac * 100, $elapsed_str, $eta_str);
        _write_progress_file($completed, $total, $frac, $elapsed_str, $eta_str);
    }
}

sub _write_progress_file {
    my ($completed, $total, $frac, $elapsed_str, $eta_str) = @_;
    if (open my $pfh, '>', "$out_prefix.progress.txt") {
        printf $pfh "Completed: %d/%d (%.1f%%)\nElapsed: %s\nETA: %s\n",
            $completed, $total, $frac * 100, $elapsed_str, $eta_str;
        close $pfh;
    }
}

sub push_recent {
    my ($name, $count) = @_;
    push @recent_taxa, { name => $name, count => $count };
    shift @recent_taxa if scalar(@recent_taxa) > $RECENT_MAX;
}

sub build_query_template {
    my $slen = "${min_slen}:${max_slen}[SLEN]";
    my $markers;
    if ($USE_CONFIG_MARKERS && @MARKER_CLAUSES) {
        $markers = scalar(@MARKER_CLAUSES) . " marker clause(s) from config";
    } else {
        $markers = "18S/28S/ITS/5.8S/H3 (built-in)";
    }
    my $content = $no_organelle
        ? $markers
        : "mito/plastid OR $markers";
    return "(SPECIES[ORGN]) AND $slen AND ($content) NOT wgs NOT tsa NOT clone NOT UNVERIFIED NOT chromosome NOT PREDICTED";
}

# ═════════════════════════════════════════════════════════════════════════════
# NCBI SUBROUTINES
# ═════════════════════════════════════════════════════════════════════════════

sub load_marker_clauses_from_regions {
    my ($path) = @_;
    return unless defined $path && -f $path;

    open my $rfh, '<', $path or do {
        warn "[NCBI-NT_Downloader] Could not open regions config [$path]: $!\n";
        return;
    };

    my $line_no = 0;
    my @clauses;
    while (my $line = <$rfh>) {
        chomp $line;
        next if $line =~ /^\s*#/;
        next unless $line =~ /\S/;

        my @cols = split /\t/, $line, 5;
        if ($line_no == 0 && defined $cols[0] && $cols[0] =~ /region_id/i) {
            $line_no++;
            next;
        }
        $line_no++;

        my ($id, $class, $enabled, $regex, $title_clause) = @cols;
        $id           //= '';
        $enabled      //= 1;
        $title_clause //= '';

        # CLI overrides take priority
        if (exists $force_enable{$id}) {
            # include regardless of enabled_default
        } elsif (exists $force_disable{$id}) {
            next;
        } else {
            next if $enabled eq '0';
        }

        next unless $title_clause =~ /\S/;
        $title_clause =~ s/^\s+|\s+$//g;
        push @clauses, "($title_clause)";
    }
    close $rfh;

    if (@clauses) {
        @MARKER_CLAUSES     = @clauses;
        $USE_CONFIG_MARKERS = 1;
    } else {
        $USE_CONFIG_MARKERS = 0;
    }
}

sub normalize_species_name {
    my ($raw) = @_;
    $raw =~ s/^\s+|\s+$//g;
    $raw =~ s/\s+/ /g;
    return $raw if $raw =~ /^txid\d+$/i;

    my @parts = split /\s+/, $raw;

    if ($name_mode eq 'as-is') {
        return $raw;
    }

    if ($name_mode eq 'trinomial') {
        return @parts >= 3 ? "$parts[0] $parts[1] $parts[2]"
             : @parts >= 2 ? "$parts[0] $parts[1]"
             : $raw;
    }

    # Legacy/default FlyGuide behavior: collapse to Genus species.
    return @parts >= 2 ? "$parts[0] $parts[1]" : $raw;
}

sub esearch_ids {
    my ($term, $email, $api_key, $retmax, $sleep, $max_retries) = @_;
    my (@ids, $count, $attempt);
    $count   = 0;
    $attempt = 0;

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
            if ($attempt < $max_retries) { sleep $sleep; next; }
            last;
        }
        last;
    }
    return ($count, @ids);
}

sub efetch_batch {
    my ($id_list, $email, $api_key, $sleep, $max_retries) = @_;
    my ($fasta, $attempt) = ('', 0);

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
        eval { $fasta = $eutil->get_Response->content; };
        if ($@) {
            if ($attempt < $max_retries) { sleep $sleep; next; }
            last;
        }
        last;
    }
    return $fasta;
}

sub format_time {
    my ($seconds) = @_;
    $seconds = int($seconds // 0);
    $seconds = 0 if $seconds < 0;
    my $h = int($seconds / 3600);
    my $m = int(($seconds % 3600) / 60);
    my $s = $seconds % 60;
    return sprintf("%02dh:%02dm:%02ds", $h, $m, $s);
}
