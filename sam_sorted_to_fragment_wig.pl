#!/usr/bin/perl
# ---------------
# sam_sorted_to_fragment_wig.pl - converts a sorted, paired-end SAM file to a variable-step WIG file of fragment
# overlap counts. Requires the read length of the raw reads used to generate the SAM file, as well as
# the minimum number of overlapping fragments you want displayed. The latter can greatly reduce the file
# size by leaving out orphaned fragments while retaining regions that are significantly enriched. Use 0
# for min_fragments to keep all reads.
# ---------------
# usage: sam_sorted_to_fragment_wig.pl sam_file wig_output read_length min_fragments
# example: ./sam_sorted_to_fragment_wig.pl csb-pgbd3_paired.sorted.sam csb-pgbd3_paired.wig 36 5
# This will create a WIG file from data generated with 36 bp reads from each end of sequenced fragments
# and will output only regions with at least 5 overlapping fragments.
# ---------------
# SAM files must be from paired-end data, and can be generated from raw reads using the read alignment program Bowtie (http://bowtie-bio.sourceforge.net/).
# SAM files must be sorted before this script is used. SAM files can be sorted using IGVTools (http://www.broadinstitute.org/igv/igvtools).
# For more information about the SAM file format, see the SAMtools page (http://samtools.sourceforge.net/).
# For a description of WIG format, see the UCSC Genome Browser description of wiggle format (https://cgwb.nci.nih.gov/goldenPath/help/wiggle.html).
# ---------------

if(@ARGV < 3) {
	
	print "usage: sam_sorted_to_fragment_wig.pl sam_file wig_output read_length min_fragments\n";
	
} else {

	$sam_in = @ARGV[0];
	$wig_out = @ARGV[1];
	$read_length = @ARGV[2] - 1;
	$min_frags = @ARGV[3];
	
	open(SAM, "<$sam_in");
	open(WIG, ">$wig_out");
	close(WIG);
	
	$last_end = 0;
	@cur_peak = ();
	$cur_chr = "chrZ";
	
	#subroutine for calculating the overlaps given a cur_peak array
	sub calc_overlaps {
		#only process if there are at least the minimum number of fragments in the peak
		if(@cur_peak >= $min_frags) {
			#the start of the peak will be the start of the first fragment in the array
			($cur_start, $z) = split(/\t/, @cur_peak[0]);
			
			#to get the end, we have to look for the largest value in the fragment ends
			$cur_end = 0;
			#cycle through all the fragments in the array
			foreach(@cur_peak) {
				#get the value for the 3' position of the fragment
				($z ,$test_end) = split(/\t/, $_);
				#if it's bigger than what's stored as $cur_end, update $cur_end with the larger value
				if($test_end > $cur_end) {
					$cur_end = $test_end;
				}
			}
			
			#get the length of the current peak
			$length = $cur_end - $cur_start;
			
			#fill an array of the peak length with 0's as a baseline
			for($x = 0; $x < $length; $x++) {
				push(@map, 0);
			}
			
			#look through each fragment and increment the map at each position the fragment overlaps
			for($j = 0; $j < @cur_peak; $j++) {
				#get the start and end of the current fragment
				($frag_start, $frag_end) = split(/\t/, @cur_peak[$j]);
				#adjust the start and end values so that the first position is 0, as in the array
				$map_start = $frag_start - $cur_start;
				$map_end = $frag_end - $cur_start;
				#for each position in the map array that the fragment overlaps,
				for($k = $map_start; $k < $map_end; $k++) {
					#increment the value of the map position
					$old_height = @map[$k];
					$map[$k] = $old_height + 1;
				}
			}
			
			open(WIG, ">>$wig_out");
			#output the overlaps at each position in the map array
			for($x = 0; $x < @map; $x++) {
				#readjust the relative position in the map to the actual genomic position
				$position = $x + $cur_start;
				#print the position and overlap to the output
				print WIG $position . "\t" . @map[$x] . "\n";
			}
			close(WIG);
			
			@map = ();
			
		}	
	}
	
	while(<SAM>) {
	
		$line = $_;
		chomp($line);
		@line = split(/\t/,$line);
		
		#Make sure this isn't a header line or from chromosome M
		if ( @line[0] eq "\@HD" || @line[0] eq "\@SQ" || @line[0] eq "\@PG" || @line[2] eq "chrM") {
		} elsif ( @line[1] == 99 || @line[1] == 163) { #if not, only process lines with flags 99 or 163.
		
			#get information about the fragment position
			$start = @line[3];
			#have to add the read length to the 3' end, because SAM format uses only the 5'-most position
			$end = @line[7] + $read_length;
			#this is the fragment line as it will be stored in the array for retrieval:
			$fragment = $start . "\t" . $end;
			
			#Check to see if this is a new chromosome
			if ( @line[2] ne $cur_chr ) {
				#if so, output the final peak from the end of the last chromosome
				
				calc_overlaps();
				
				#and start processing the new chromosome.
				$cur_chr = @line[2];
				print "Now processing " . $cur_chr . "\n";
				
				#output the header line for the new chromosome to the wig file
				open(WIG, ">>$wig_out");
				print WIG "variableStep chrom=" . $cur_chr . "\n";
				close(WIG);
				
				#reset the end point position and current peak array
				$last_end = 0;
				@cur_peak = ();
			
			} else {			
				#if this fragment is outside the last peak, output the last peak and start storing
				#fragments for a new peak.
				if($start > $last_end) {
				
					calc_overlaps();
				
					@cur_peak = ();
				}
				
			}
			
			push(@cur_peak, $fragment);
			
			if($end > $last_end) {
				$last_end = $end;
			}
		}
	}
	
	#output the last peak
	
	calc_overlaps();
	
	close(SAM);

}