
# RNA-seq_noise_reduction.pl version 1

# This script was created by Andre E Minoche

# last upadte April 26, 2015

# Please cite this Genome Biology article



######## Description of the algorithm
# The default settings reduce the overall RNA-seq coverage by 10% [a] of the local peak coverage
# (95 percentile [b], in 1-kbp windows [c]). If the coverage difference between two intron hints
# (gff format) is greater than 90% the intron hint with the lower coverage gets removed. The
# mRNA-seq coverage (wig format) gets reduced within boundaries of introns by 50% [d] of the
# adjacent exon coverage. Introns are considered if they are smaller or equal than 50 kbp in
# size and show a coverage drop of at least 50% at their exon-intron junction when comparing
# the coverage 10 bases upstream and downstream of each junction.


our $ExCutOff=0.10; # [a]
our $ExMaxVal=0.95; # [b]
our $ExonFlanking=500; # [c]
our $cleanIntronBy=0.5; # [d]

# [c] also triggers the update frequency for which the local coverage peak is calculated
our $updateMaxExonCalcEveryXBases=$ExonFlanking*0.1;

our @Ap2v, @Ap, %IntrStEn, %Int2Mult, $cScaf;

print STDERR "read intron hints...\n";
open(IN1, "<$ARGV[0]") || die "can not read $ARGV[0]";
while(<IN1>){chomp;@_=split("\t",$_);

  if($_[8] =~/^mult=([0-9]+)/){$mult=$1}else{next;$mult=1}
  $IntrStEn{$_[0]}{$_[3]}{$_[4]}=$_;
  $Int2Mult{$_[0].":".$_[3].":".$_[4]}=$mult;

}close(IN1);

# -----------------------------------------------
# debug
# to see if it prints the content

# foreach $x (keys %IntrStEn) {
#   foreach $y (keys $IntrStEn{$x}) {
#     foreach $key (keys $IntrStEn{$x}{$y}) {
#       print "$IntrStEn{$x}{$y}{$key}\n";
#     }
#   }
# }
# -----------------------------------------------

open(INTR, ">$ARGV[2]") || die "can not write to $ARGV[2]";
open(EXON, ">$ARGV[3]") || die "can not write to $ARGV[3]";

print STDERR "read wig file containing the exons coverage...\n";
open(IN1, "<$ARGV[1]") || die "can not read $ARGV[1]";
while(<IN1>){

  chomp;
  if(/^[^0-9]/){

  	if(/chrom/){
        # previous version by Andre Minoche
        # the split is performed on the "="
        # however, some studies produce chromosome names that carry the length
        # ... separated by a "=" to the chromosome name
        # ... such as:
        # Bvchr5.seq2_len=11159305
        #
  	    # @t=split("=",$_);
        #
        # new version (Matteo Schiavinato)
        # this new version splits on the "chrom=" pattern
        # ... to make sure that the whole chromosome name that follows
        # ... is retained
        #
        @t=split("chrom=",$_);
        #
  	    processScaf() if @Ap>0;
  	    undef(@Ap2v);
  	    undef(@Ap);
  	    $cP=0;
  	    $cScaf=$t[1];
  	 }
  	 print EXON $_."\n";
  	 next;
  }

  ($p,$v)=split(" ",$_);

  push @Ap,$p; # array that contains the positions included in the wig file for the current scaffold

  # create an array that represent the coverage of the current scaffold
  # positions not included in the wig file have coverage = 0
  while($cP<=$p){$Ap2v[$cP]+=0; $cP++;}
  $Ap2v[$p]+=$v;

}close(IN1);

processScaf() if @Ap>0;


sub processScaf{

  # create copies of the coverage representation, that will be modified
  @Ap2vC=@Ap2v;
  @Ap2vI=@Ap2v;

  $lSt=0;$lEn=0;
  #print STDERR "remove low coverage intron hints for $cScaf ...          \r";
  print STDERR "remove low coverage intron hints for $cScaf ...          \n";

  $del=1;$countDelIntrons=0;
  while($del==1){ # repeat the comparison of adjacent intron hints until no further intron hints got deleted
	  $del=0;
	  # for each adjacent intron, test if they overlap, and if the coverage difference meets the threshold
	  foreach $cSt (sort {$a <=> $b} keys %{$IntrStEn{$cScaf}}){
		foreach $cEn (sort {$a <=> $b} keys %{$IntrStEn{$cScaf}{$cSt}}){

		  if (overlap($cSt,$cEn,$lSt,$lEn)>0){
			$cMult=$Int2Mult{$cScaf.":".$cSt.":".$cEn};
			$lMult=$Int2Mult{$cScaf.":".$lSt.":".$lEn};

			if($cMult*0.1>$lMult){
			  $del=1; $countDelIntrons++;
		# 	  print STDERR "del intron: $cScaf $lSt $lEn, ".($cMult*0.1)."> $lMult\n  ".$IntrStEn{$cScaf}{$cSt}{$cEn}."\n  ".$IntrStEn{$cScaf}{$lSt}{$lEn}."\n";
			  delete($IntrStEn{$cScaf}{$lSt}{$lEn});
			}elsif($lMult*0.1>$cMult){
			  $del=1; $countDelIntrons++;
		# 	  print STDERR "del intron: $cScaf $cSt $cEn, $cMult <".($lMult*0.1)."\n  ".$IntrStEn{$cScaf}{$cSt}{$cEn}."\n  ".$IntrStEn{$cScaf}{$lSt}{$lEn}."\n";
			  delete($IntrStEn{$cScaf}{$cSt}{$cEn});
			  next;
			}
		  }
		  $lSt=$cSt; $lEn=$cEn;
		}
	  }
  }
  #print STDERR "deleted intron hints $countDelIntrons         \r";
  print STDERR "deleted intron hints $countDelIntrons         \n";

  # reduce RNA-seq coverage in intron regions using cleaned intron hints.
  #print STDERR "clean intron regions for $cScaf...                  \r";
  print STDERR "clean intron regions for $cScaf...                  \n";
  foreach $cSt (sort {$a <=> $b} keys %{$IntrStEn{$cScaf}}){
    foreach $cEn (sort {$a <=> $b} keys %{$IntrStEn{$cScaf}{$cSt}}){

		$mult=$Int2Mult{$cScaf.":".$cSt.":".$cEn};
		$ExonCov=($Ap2v[($cSt-10)]+$Ap2v[($cEn+10)])/2;

		$mult=$Int2Mult{$cScaf.":".$cSt.":".$cEn};
		for $ip ($cSt..$cEn){ $Ap2vI[$ip]+=$mult; } # increase the Wig coverage by the number of split reads representing an intron hint. For later when the exon coverage gets adjusted.
		print INTR $IntrStEn{$cScaf}{$cSt}{$cEn}."\n"; # print the cleaned intron hint

		# clean the unspliced introns if intron < 50,000
		# and if the coverage drop is at least 50% compared to each adjacent exon borders
		if ($Ap2v[($cSt-10)]/($Ap2v[($cSt+10)]+0.001)<2 or $Ap2v[($cEn+10)]/($Ap2v[($cEn-10)]+0.001)<2){
	# 	  print STDERR "no cleaning: $cScaf $cSt $cEn, ".$Ap2v[($cSt-10)]."-".($Ap2v[($cSt+10)]+0.001)." vs ".$Ap2v[($cEn+10)]."-".($Ap2v[($cEn-10)]+0.001)." $mult  \n";
		  next;
		}else{
	# 	  print STDERR "----cleaning: $cScaf $cSt $cEn, ".$Ap2v[($cSt-10)]."-".($Ap2v[($cSt+10)]+0.001)." vs ".$Ap2v[($cEn+10)]."-".($Ap2v[($cEn-10)]+0.001)." $mult \n";
		}

		next if abs($cEn-$cSt+1)>50000;

		$reduceCovBy=sprintf("%.0f",($ExonCov*$cleanIntronBy));
	# 	print STDERR "rm $cSt - $cEn\n";
		for($cSt..$cEn){

			  $Ap2vC[$_]=($Ap2vC[$_]>$reduceCovBy)? ($Ap2vC[$_]+($reduceCovBy*(-1))):0;
		}
    }
  }

  # reduce the RNA-seq coverage in wig, by 10% of the local maximum (95 percentile, ExMaxVal)
  $positionCount=0;
  #print STDERR "adjust exon coverage for $cScaf...            \r";
  print STDERR "adjust exon coverage for $cScaf...            \n";
  for ($i=0; $i<=$#Ap; $i++){ # for each position int he wig file
    $cVal=$Ap2vC[$Ap[$i]];
    next if !$cVal>0;
    $cStmp="";

    # test if previous position in wig file is == curent position -1
    $cPosIsAdjecent=0;
    if (($Ap[$i]-1)==$Ap[($i-1)]){
      $positionCount++;
      $cPosIsAdjecent=1;
    }else{
      $positionCount=0;
    }

    # adjust the local exon coverage threshold if current base is the first base after gap or every $updateMaxExonCalcEveryXBases bases, or first base in scaffold
    if($cPosIsAdjecent==0 or $positionCount%$updateMaxExonCalcEveryXBases==1 or $i==0){
      $min=(($Ap[$i]-$ExonFlanking)<1)? 1:($Ap[$i]-$ExonFlanking);
      $max=(($Ap[$i]+$ExonFlanking)>$Ap[$#Ap])? $Ap[$#Ap]:($Ap[$i]+$ExonFlanking);
      @tmp=sort {$a <=> $b} grep {$_%10==1} @Ap2vI[$min..$max];
      $cMaxVal=$tmp[(sprintf("%.0f",(scalar(@tmp)*$ExMaxVal))-1)];
      $cStmp.="r ";
      $reduceCovBy=sprintf("%.0f",($ExCutOff*$cMaxVal));
    }

    $cVal=($cVal>$reduceCovBy)? ($cVal+($reduceCovBy*(-1))):0; # reduce the current coverage value, by reduceCovBy, but not less than 0

    $cStmp.="min: $min, max: $max, cMax: $cMaxVal, array size ".scalar(@tmp)." reduceCovBy: $reduceCovBy";

    # only print positions with a coverage > 0
    if ($cVal==0){
# 		print EXON "$Ap[$i]\t".$cVal."\t".$Ap2v[$Ap[$i]]."* $cStmp\n";
		next;
    }else{
# 		print EXON "$Ap[$i]\t".$cVal."\t".$Ap2v[$Ap[$i]]." $cStmp\n";
		print EXON "$Ap[$i]\t".$cVal."\n";
    }
  }

}
close(INTR);
close(EXON);



# return the length two fragments are overlapping
sub overlap{
	my($R1,$R2,$Q1,$Q2)=@_;
	($QuerB,$QuerE) = sort {$a <=> $b} ($Q1,$Q2);
	($RefB,$RefE)   = sort {$a <=> $b} ($R1,$R2);
	$ovlLen=(sort {$a <=> $b} ($QuerE,$RefE))[0]-(sort {$b <=> $a} ($QuerB,$RefB))[0]+1; # ovl if positive
	$returnVal=($ovlLen>0)? $ovlLen:0;
	return $returnVal;
}
