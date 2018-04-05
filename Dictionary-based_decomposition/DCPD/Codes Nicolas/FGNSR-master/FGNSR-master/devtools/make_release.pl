#!/usr/bin/perl

use warnings;
use strict;
use File::Copy "cp";
use File::stat;
use Cwd;

my $version;
if (scalar @ARGV != 1) {
    print "Usage: make_release.pl <version>\n";
    exit(0);
}
else {
    $version = $ARGV[0];
}

# Files are taken from pwd
my $co_dir = getcwd();

# Files are put into this directory
my $target_dir = '.';

my $package_dir = $target_dir . '/fgnsr-' . $version;

if (-e $package_dir) {
    die "Error: $package_dir already exists.\nRemove old directory first.\n";
}

# Define all directories in the package.
my @new_directories = (
    "$package_dir",
    "$package_dir/matlab",
    "$package_dir/matlab/private",
    "$package_dir/ckernel",
);
print "Creating directory tree for $package_dir ...\n";
for my $dir (@new_directories) {
    my $success = mkdir($dir);
    if (not $success) {
        die "Could not create directory $dir because $!\n";
    }
}

my %filemap = (
    'LICENCES'                       => 'dist/LICENSES.dist',
    'README.md'                      => 'README.md',
    'CHANGES'                        => 'dist/CHANGES.dist',
    'setup_fgnsr.m'                  => 'setup_fgnsr.m',
    'ckernel/fgnsr_ss.h'             => 'ckernel/fgnsr_ss.h',
    'ckernel/fgnsrdef.h'             => 'ckernel/fgnsrdef.h',
    'ckernel/util.h'                 => 'ckernel/util.h',
    'ckernel/util.c'                 => 'ckernel/util.c',
    'ckernel/proj.c'                 => 'ckernel/proj.c',
    'ckernel/proj.h'                 => 'ckernel/proj.h',
    'ckernel/parproj_priv.c'         => 'ckernel/parproj_priv.c',
    'matlab/parproj.m'               => 'matlab/parproj.m',
    'matlab/fgnsr.m'                 => 'matlab/fgnsr.m',
    'matlab/example_midpoints.m'     => 'matlab/example_midpoints.m',
    'matlab/proj_omega.m'            => 'matlab/proj_omega.m',
    'matlab/private/fgnsr_init.m'    => 'matlab/private/fgnsr_init.m',
    'matlab/private/FastSepNMF.m'    => 'matlab/private/FastSepNMF.m',
    'matlab/private/nnlsHALSupdt.m'  => 'matlab/private/nnlsHALSupdt.m',
);
print "Copying files from $co_dir to $package_dir\n";
for my $target (keys %filemap) {
    my $targetfile = "$package_dir/$target";
    my $srcfile    = "$co_dir/$filemap{$target}";

    my $success =  cp($srcfile, $targetfile);
    if (not $success) {
        die "Error: Could not copy $srcfile to $targetfile because $!\n";
    }

    $success = chmod((stat($srcfile)->mode), $targetfile);
    if (not $success) {
        warn "Warning:  could not set permissions for $targetfile\n";
    }
}
