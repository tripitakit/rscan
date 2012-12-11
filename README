rscan - version: 0.1
=====================

Performs nucleotide scoring and color masking in multiple sequence alignments, to help in DNA signature search.
Runs within irb, in a terminal emulator with ANSI color support.
Requires gems: bio, rainbow.

Usage:
------
  ./rscan

Input:
------
Multiple DNA-sequence alignment fasta-formatted file. Groups definitions.
Default scoring formula =  1 - (ka/5)*(1-a) - (kb/10)*b
	 a: frequency of nt_base in the ingroup (consenus)
	 b: frequency of nt_base in the outgroup  (aspecificity)
	 ka: consensus coefficient : 20 :strict, 10 :high (majority rule high), 3 :low (majority rule low)
	 kb: aspecificity tollerance coefficient : 0 :allow , 10 :penalty (allowed with penality), 20 :forbid (not allowed)


IRB user interface command list:
=================================

Alignment and groups
--------------------
  open "filename"
     open a fasta input file and create the working alignment
     usage:
     >> open "test.fasta"
  
  labels?
     print the sequences' index and label
     usage:
     >> labels?
    
  groups?
     print groups array
     usage:
     >> groups?

  set_groups [arrays and/or ranges]
     set the groups of sequences, expects an array of arrays or ranges
     usage, with a 16 sequences alignemnt:
     >> set_groups [0..4, 5..9, 10..15]
     >> set_groups [0..4, [5,6,7], 8..14, [15]]
  
  
Set Consensus coefficient (ka)
------------------------------
  con :preset_value
     use a preset value for consensus coefficient ka
     expects a symbol [:strict | :high | :low]
     usage:
     >> con :strict 
    
  ka numeric_value
     setter for ka (consensus coefficient) 
     expects a numeric value 
     (standard values for default formula are: 20, 10, 3)
     usage:
     >> ka 10
    
    
Set Aspecificity-tollerance coefficient (kb)
-------------------------------------------- 
  asp :preset_value
     use a preset value for aspecificity tollerance coefficient kb
     expects a symbol [:forbid | :penalty | :allow]
     usage:
     >> asp :allow

  kb numeric_value
     setter for kb (aspecificity tollerance coefficient)
     expects an integer value
     (standard values for default formula are: 20, 10, 0)
     usage:
     >> kb 10


Scan and Scores
---------------
  scan
     scan the alignment with the current parameters and print the score-masking color shading
     usage:
     >> scan
       
  color_ranges [n1, n2, n3, n4]
     defines the borders values of score ranges, for color masking
     expects an array of 4 values, default is [0.5, 0.7, 0.8, 0.9]
     ranges are defined as:
       score < n1
       n1 <= score < n2
       n2 <= score < n3
       n4 <= score < n4
       score >= n4
     usage:
     >> color_ranges [0.3, 0.45, 0.6, 0.88]
    
  set_formula "a function of a, b, @ka, @kb"
     setter for the scoring formula
     expects a string; variables and coefficient: a, @ka; b, @kb
     usage:
     >> setformula "1 - (a * @ka) - (b * @kb)"  

  formula?
     print the scoring formula
     usage:
     >> formula?
    

Utilities
---------
  clear
     clear the console screen
     usage:
     >> clear
    
  test
     perform a test scan 
     usage:
     >> test	
