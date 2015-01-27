#!/usr/bin/perl
# M. Toups
# 12/1/14
#
# Short script to submit optical library generation with
# project.py.
#
use strict;
use warnings;
use File::Basename qw(fileparse);
use Cwd qw(abs_path);

my ($xml, $fcl, $workdir, $check, $merge) = @ARGV;

if(!(defined($xml) || defined($fcl) || defined($workdir))) {
    print "Usage: perl SubmitCommand.pl <project.py xml file> <fcl file> <work dir> [<checkana?>] [<merge?>]\n";
    exit(1);
}

if(!(-e "${workdir}/buildopticallibrary")) {
    mkdir "${workdir}/buildopticallibrary", 0775;
}

if(!(-e "${workdir}/buildopticallibrary/xml")) {
    mkdir "${workdir}/buildopticallibrary/xml", 0775;
}

if(!(-e "${workdir}/buildopticallibrary)/fcl")) {
    mkdir "${workdir}/buildopticallibrary/fcl", 0775;
}

my $NPhotonsPerVoxel = 30000;
my $NVoxelsPerJob = 240;
my $NJobs = 9375;
my $StartJob = 0;

for(my $i=$StartJob; $i<$StartJob+$NJobs; $i++) {

    my ($fcl_filename, $fcl_directories, $fcl_suffix) = fileparse($fcl,".fcl");
    my ($xml_filename, $xml_directories, $xml_suffix) = fileparse($xml,".xml");

    if(!(-e "${workdir}/buildopticallibrary/fcl/${fcl_filename}_${i}${fcl_suffix}")) { 

	open IN, "${fcl}" or die $!;
	open OUT, ">${workdir}/buildopticallibrary/fcl/${fcl_filename}_${i}${fcl_suffix}" or die $!;
	
	while(<IN>) {
	    if(/physics.producers.generator.FirstVoxel/) {
		print OUT "physics.producers.generator.FirstVoxel: " . $NVoxelsPerJob*$i . "\n";
	    }
	    elsif(/physics.producers.generator.LastVoxel/) {
		print OUT "physics.producers.generator.LastVoxel: " . ($NVoxelsPerJob*($i+1)-1) . "\n";
	    }
	    elsif(/physics.producers.generator.N/) {
		print OUT "physics.producers.generator.N: ${NPhotonsPerVoxel}\n";
	    }
	    else {
		print OUT "$_";
	    }
	}
	close IN or die $!;
	close OUT or die $!;
    }

    if(!(-e ">${workdir}/buildopticallibrary/xml/${xml_filename}_${i}${xml_suffix}")) {

	open IN, "${xml}" or die $!;
	open OUT, ">${workdir}/buildopticallibrary/xml/${xml_filename}_${i}${xml_suffix}" or die $!;

	while(<IN>) {
	    if(/ENTITY name/) {
		print OUT '<!ENTITY name "gen_photon_ball_' . $i . '">' . "\n";
	    }
	    elsif(/\<fcl\>(.*)\<\/fcl\>/) {
		print OUT "    \<fcl\>" . abs_path(${workdir}) . "\/buildopticallibrary\/fcl\/${fcl_filename}_${i}${fcl_suffix}\<\/fcl\>\n";
	    }
	    else {
		print OUT "$_";
	    }
	}

	close IN or die $!;
	close OUT or die $!;
    }

    if(!defined($check)) {
	if(!defined($merge)) {
	    system qq|project.py --xml ${workdir}/buildopticallibrary/xml/${xml_filename}_${i}${xml_suffix} --submit|;
	    sleep 1;
	} else {
	    system qq|project.py --xml ${workdir}/buildopticallibrary/xml/${xml_filename}_${i}${xml_suffix} --merge|;
	}
    } else {
	system qq|project.py --xml ${workdir}/buildopticallibrary/xml/${xml_filename}_${i}${xml_suffix} --checkana|;
    } 

}


