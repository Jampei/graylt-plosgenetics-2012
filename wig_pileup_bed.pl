#!/usr/bin/perl
# wig_pileup.pl - Reads a wiggle file to get the cumulative fragment pileup of reads
# over a set of regions defined by a start site, a direction, and the number of 
# bases up and downstream from the start sites that should be analyzed.
# usage: wig_pileup_bed.pl in_locations in_wig upstream downstream out_file

if (@ARGV != 5) {
	print "usage: wig_pileup_bed.pl in_locations in_wig upstream downstream out_file \n";
} else {

	#declare input variables
	$in_loc = @ARGV[0];
	$in_wig = @ARGV[1];
	$up = @ARGV[2];
	$down = @ARGV[3];
	$out_file = @ARGV[4];
	
	#Generate an array to store the pileup values in, starting out at zero for each position
	$total_region = $up + $down;
	
	for($i = 0; $i < $total_region; $i++) {
		push(@pileup, 0);
	}
	
	#calculate the start and end positions for each location to make searching through the 
	#wiggle file simpler later
	
	open(LOC, "<$in_loc");
	while(<LOC>) {
		
		$line = $_;
		chomp($line);
		
		@line_split = split(/\t/, $line);
		
		if(@line_split[5] eq "+") {
			$start = @line_split[1] - $up;
			$end = @line_split[1] + $down;
		} else {
			$start = @line_split[2] - $down;
			$end = @line_split[2] + $up;
		}
		
		$new_line = $start . "\t" . $end . "\t" . @line_split[5];
		
		#put the adjusted locations in arrays for each chromosome
		push(@{"@line_split[0]"}, $new_line);

	}
	close(LOC);
	
	#then parse through the wiggle file to check if each line falls within a region of interest
	$display_counter = 0;
	$total_counter = 0;
	$cur_chr = "Locations";
	open(WIG, "<$in_wig");
	while(<WIG>) {
		
		$line = $_;
		chomp($line);
		
		#check to see if the line is a chromosome header
		@chr_check = split(/=/,$line);
		
		if (@chr_check[0] eq "variableStep chrom") {
			print "Finished processing $cur_chr.\n";
			$cur_chr = @chr_check[1];
			
		} else {
			
			#if it's not a header, check through the regions for the current chromosome
			
			for ($i = 0; $i < @{"$cur_chr"}; $i++) {
				
				@region_split = split(/\t/, @{"$cur_chr"}[$i]);
				
				#check to see if the wiggle location falls in the region
					
				@split_wig = split(/\t/, $line);
				
				if ( @split_wig[0] < @region_split[0] ) {
					#if the wiggle position is below the start of the current region
					#stop looping through regions. This should improve speed.
					last;
				} elsif (@split_wig[0] > @region_split[1]) {
					#if the wiggle position is above the end of the current region, the region
					#should be removed from the array to speed comparisons
					shift(@{"$cur_chr"});
					$i--;
				} elsif ( @split_wig[0] <= @region_split[1] ) {
					#if it does fall in the region, figure out what position it's in
					#this will be dependent on the direction of the region (@region_split[3])
					
					if (@region_split[2] eq "+") {
						$corrected_position = @split_wig[0] - @region_split[0];
					} elsif (@region_split[2] eq "-") {
						$corrected_position = @region_split[1] - @split_wig[0];
					}
					
					#finally, increase the value of the corrected position by the
					#height value for the position from the wiggle file.
					$new_height = @pileup[$corrected_position] + @split_wig[1];
					splice(@pileup,$corrected_position,1,$new_height);
					
				}
				
			}
			
		}
		
		$display_counter++;
		$total_counter++;
		
		if ($display_counter == 10000000) {
			print "Processed $total_counter wiggle lines\n";
			$display_counter = 0
		}
		
	}
	
	#After going through the wiggle file, we should be able to output the overlaps
	
	open(OUTPUT, ">$out_file");
	
	for ($i = 0; $i < @pileup; $i++) {
		
		$position = 0 - $up + $i;
		
		print OUTPUT $position . "\t" . @pileup[$i] . "\n";
		
	}
	
}

