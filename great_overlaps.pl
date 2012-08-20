#!/usr/bin/perl
#great_overlaps.pl - calculates the number of overlaps between ChIP-seq peaks and GREAT regulatory
#domain files for use with the GREAT calculateBinomialP program. The position used for the peak is
#the center of the region provided in the BED file.
#usage: great_overlaps.pl domains.rd peaks.bed

#set variables for input and output
$rd_file = @ARGV[0];
$summ_file = @ARGV[1];

open(RD, "<$rd_file");

#load the regulatory domains for comparison with the peaks
while(<RD>) {

	$line = $_;
	chomp($line);
	
	@line_array = split(/\t/, $line);
	
	$rd_loc = @line_array[1] . "\t" . @line_array[2];
	
	push(@{"@line_array[0]"},$rd_loc);

}

close(RD);
#open the peak file and check each peak to see if it hits the regulatory domains.

open(PEAKS, "<$summ_file");

$hits = 0;
$total_peaks = 0;

while(<PEAKS>) {

	$line = $_;
	chomp($line);

	@peak_split = split(/\t/, $line);
	
	$peak_length = @peak_split[2] - @peak_split[1];
	$summit_relative = $peak_length / 2;
	$summit_position = @peak_split[1] + $summit_relative;
	
	$summ_hit = 0;

	for($i=0; $i < @{"@peak_split[0]"}; $i++) {
		
		@rd_loc = split(/\t/,@{"@peak_split[0]"}[$i]);
		
		#first, check to see if the chromosome matches
		
		if(@rd_loc[0] > $summit_position) {
			last;
		} elsif (@rd_loc[1] < $summit_position - 50000) {
			shift(@{"@peak_split[0]"});
			$i--;
		} elsif ($summit_position >= @rd_loc[0] && $summit_position <= @rd_loc[1]) {
			
			#if that's true, and this peak hasn't hit a region yet 
			#then we numerate the hits variable
			
			if($summ_hit == 0) {
			
				$hits++;
				
				$summ_hit = 1;
				
			}
		}
	}
	
	#we need to keep track of the total number of peaks, too
	$total_peaks++;

}

close(PEAKS);

#print the total peaks and the number of hits
print $total_peaks . "\t" . $hits;