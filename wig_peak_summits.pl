#!/usr/bin/perl
# wig_peak_summits.pl - looks through a wiggle file for the highest points in a given set of bed regions.
# These are the areas of highest fragment overlap - height, start, end, length, and center are
# calculated and output.
# usage: wig_peak_summits.pl in_wig in_peaks out_summits

if (@ARGV != 3) {
	print "usage: wig_peak_summits.pl in_wig in_peaks out_summits\n";
} else {
	
	$wig_in = @ARGV[0];
	$peaks_in = @ARGV[1];
	$out_file = @ARGV[2];
	
	open(PEAKS, "<$peaks_in");
	while(<PEAKS>) {
		
		$line = $_;
		chomp($line);
		
		@line_split = split(/\t/, $line);
		
		push(@{"@line_split[0]"}, $line);

	}
	close(PEAKS);
	
	sub process_chr {
		print "Finished reading $cur_chr.\n";
			#if it is, check through the regions for that chromosome to do the calculations.
			
			for ($i = 0; $i < @{"$cur_chr"}; $i++) {
				
				#split up the region. 0 is chr, 1 is start, 2 is end.
				@region_split = split(/\t/, @{"$cur_chr"}[$i]);
				
				#set the variables to keep track of
				$high_start = 0;
				$high_end = 0;
				$high_value = 0;
				$find_end = 0;
				
				#look through each wiggle point for the chromosome
				for ($j = 0; $j < @chr_wig; $j++) {
					
					#split the wiggle line. 0 is position, 1 is height
					@split_wig = split(/\t/, @chr_wig[$j]);
					
					if ( @split_wig[0] < @region_split[1] ) {
						#if the wiggle position is below the start of the current region
						#get rid of the wiggle position so the search is faster
						shift(@chr_wig);
						$j--;
					} elsif (@split_wig[0] > @region_split[2]) {
						#if the wiggle position is above the end of the current region, 
						#stop looking through the wiggle data
						
						last;
						
					} elsif ( @split_wig[0] <= @region_split[2] ) {
						#if it does fall in the region, check for top position
						
						if (@split_wig[1] > $high_value) {
							#if this is a new high point, set the location, the height, and
							#set the variable to check for where it ends
							$high_start = @split_wig[0];
							$high_value = @split_wig[1];
							$find_end = 1;
							
						} elsif ( ($find_end == 1) && (@split_wig[1] < $high_value) ) {
							#if a new high point has been found, and the current location is 
							#lower than the high point, it's the end of the overlap plateau.
							$high_end = @split_wig[0] - 1;
							$find_end = 0;
							
						}
						
					}
				
				}
			
			#in some cases, the highest overlap for the peak is at the end of the region.
			#when that happens, the end won't be found properly, and should be the end of the region.
			
			if ($high_end < $high_start) {
				$high_end = @region_split[2];
			}
			
			#after looking through the wiggle positions, print the top point for that peak.

			$summit_length = $high_end - $high_start + 1;
			$center = int($high_start + $summit_length / 2);
			$rel_center = $center - @region_split[1];
			$peak_length = @region_split[2] - $region_split[1];
			$rel_center_ratio = $rel_center / $peak_length;
			
			print OUTPUT @{"$cur_chr"}[$i] . "\t" . $high_start . "\t" . $high_end . "\t" . $high_value . "\t" . $summit_length . "\t" . $center . "\t" . $rel_center . "\t" . $rel_center_ratio . "\n";
			
			}
			
			@chr_wig = ();
			$cur_chr = @chr_check[1];
			
	}
	
	#Then, look through the wiggle file for positions that are in the regions
	$display_counter = 0;
	$total_counter = 0;
	
	$cur_chr = "Regions";
	
	open(WIG, "<$wig_in");
	open(OUTPUT, ">$out_file");

	#Read in the wiggle file. I'll do this a chromosome at a time, then process that chromosome.
	while(<WIG>) {
		
		$line = $_;
		chomp($line);
		
		#check to see if the line is a chromosome header
		@chr_check = split(/=/,$line);
		
		if (@chr_check[0] eq "variableStep chrom") {
			#if this is a new chromosome, run the process chromosome subroutine, above.
			process_chr();
			
		} else {
			
			push(@chr_wig, $line);			
			
		}
		
		$display_counter++;
		$total_counter++;
		
		if ($display_counter == 10000000) {
			print "Processed $total_counter wiggle lines\n";
			$display_counter = 0;
		}
		
	}
	
	#we'll need to run the process_chr one more time in order to get the last chromosome output:
	process_chr();
	
	close(WIG);
	close(OUTPUT);
	
	
}
	

