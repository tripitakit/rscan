# rscan.rb
# 
# usage: irb --simple-prompt -r rscan.rb
# require a terminal emulator whit ansi color support
# tested on: OSX terminal.app, GNOME-terminal
# author email: patrick.demarta@gmail.com

require "./signs"

def clear
  system "clear"
  puts "rscan 0.1 :: Alignment shader for signature-sequence search.".color(:blue)
  puts "\n[QuickHelp: type man for the list of commands]".color(:magenta)
  puts "\n"
  return Time.now
end

def open(file)
  # open a fasta input file, create the working alignment
  # usage:
  # >> open "test.fasta"
  @ali = Signs.new(file)
  return @ali
end

def set_groups(an_array)
  # set the groups of sequences, expects an array of arrays or ranges
  # usage:
  # >> set_groups [0..4, 5..9, 10..15]
  # >> set_groups [0..4, [5,6,7], 8..14, [15]]
  @ali.groups = an_array
  return :done
end

def groups?
  # print groups array
  # usage:
  # >> groups
  @ali.groups
  return :done
end

def labels?
  # print the sequences labels
  # usage:
  # >> sequences
  @ali.labels
  return :done
end

def scan
  # scan the alignment with the current parameters and print the score-masking color shading
  # usage:
  # >> scan
  @ali.scan
  @ali.print_scores
  return :done
end

def ka(int)
  # setter for ka (consensus coefficient) 
  # expects a numeric value
  # usage:
  # >> ka 10
  @ali.ka = int
  return :done
end

def kb(int)
  # setter for kb (aspecificity tollerance coefficient)
  # expects a numeric value
  # usage:
  # >> kb 10
  @ali.kb = int
  return :done
  
end

def con(preset)
  # presets for consensus coefficient
  # expects a symbol [:strict | :high | :low]
  # usage:
  # >> con :strict 
  case preset
  when :strict then @ali.ka = 20
  when :high then @ali.ka = 10
  when :low then @ali.ka = 3
  end
end  

def asp(preset)
  # presets for aspecificity tollerance coefficient
  # expects a symbol [:forbid | :penalty | :allow]
  # usage:
  # >> asp :allow
  case preset
  when :forbid then @ali.kb = 20
  when :penalty then @ali.kb = 10
  when :allow then @ali.kb = 0
  end
end    


def color_ranges(*ary)
  # defines the borders values of score ranges, for color masking
  # expects a list of 4 values
  if ary.size == 4
    @ali.score_color_ranges = ary
  else
    puts "4-values list is required, in example:"
    puts ">> color_ranges 0.2, 0.4, 0.6, 0.8)"
    @ali.score_color_ranges =  [0.5, 0.7, 0.8, 0.9]
    puts "Reset to default: 0.5, 0.7, 0.8, 0.9"
  end  
    return :done
  end

def set_formula(formula)
  # setter for the scoring formula
  # expects a string; variables and coefficient: a, ka; b, kb
  # usage:
  # >> set_formula "1 - (a * ka) - (b * kb)"
  @ali.scoring_formula = formula.gsub("k", "@k") # set ka and kb as instance variables
  return :done
end

def formula?
  # print the scoring formula
   @ali.scoring_formula.gsub("@k", "k") # hide the instance variables symbol "@" in @ka and @kb
end

def test
  consenus = [:strict, :low]
  aspecificity = [:forbid, :penalty]
  open "test.fasta"
  set_groups [0..7,8..12,13..20,21..26]
  groups?
  params = consenus.product aspecificity
  params.each do |pair|
    con pair[0]
    asp pair[1]
    print "Consensus: #{pair[0]}(#{@ali.ka}), Aspecificity: #{pair[1]}(#{@ali.kb})}\n"
    scan 
    print "\n"
  end
  return :done
end

def man
  manual = <<EOF 
  
rscan - version: 0.1

  Performs nucleotide scoring and color masking in multiple sequence alignments, to help in DNA signature search.
  Runs within irb, in a terminal emulator with ANSI color support.
  Requires gems: bio, rainbow.


  Usage:  irb -r rscan


  Required input: a fasta formatted  multiple sequence (DNA) alignment. Groups definitions.

  Default scoring formula =  1 - (ka/5)*(1-a) - (kb/10)*b
  	* a: frequency of nt_base in the ingroup (consenus)
  	* b: frequency of nt_base in the outgroup  (aspecificity)
  	* ka: consensus coefficient : 20 :strict, 10 :high (majority rule high), 3 :low (majority rule low)
  	* kb: aspecificity tollerance coefficient : 0 :allow , 10 :penalty (allowed with penality), 20 :forbid (not allowed)


  IRB user interface command list:
  =================================

  Alignment and groups
  --------------------
    open "filename"
      # open a fasta input file and create the working alignment
      # usage:
      # >> open "test.fasta"

    labels?
      # print the sequences' index and label
      # usage:
      # >> labels?

    groups?
      # print groups array
      # usage:
      # >> groups?

    set_groups [arrays and/or ranges]
      # set the groups of sequences, expects an array of arrays or ranges of sequence's indexes
      # usage, with a 16 sequences alignemnt:
      # >> set_groups [0..4, 5..9, 10..15]
      # >> set_groups [0..4, [5,6,7], 8..14, [15]]


  Set Consensus coefficient (ka)
  ------------------------------
    con :preset_value
      # use a preset value for consensus coefficient ka
      # expects a symbol [:strict | :high | :low]
      # usage:
      # >> con :strict 

    ka numeric_value
      # setter for ka (consensus coefficient) 
      # expects a numeric value 
      # (standard values for default formula are: 20, 10, 3)
      # usage:
      # >> ka 10


  Set Aspecificity-tollerance coefficient (kb)
  -------------------------------------------- 
    asp :preset_value
      # use a preset value for aspecificity tollerance coefficient kb
      # expects a symbol [:forbid | :penalty | :allow]
      # usage:
      # >> asp :allow

    kb numeric_value
      # setter for kb (aspecificity tollerance coefficient)
      # expects an integer value
      # (standard values for default formula are: 20, 10, 0)
      # usage:
      # >> kb 10


  Scan and Scores
  ---------------
    scan
      # scan the alignment with the current parameters and print the score-masking color shading
      # usage:
      # >> scan

    color_ranges [n1, n2, n3, n4]
      # defines the borders values of score ranges, for color masking
      # expects an array of 4 values, default is [0.5, 0.7, 0.8, 0.9]

      # ranges are defined as:
      #   score < n1
      #   n1 <= score < n2
      #   n2 <= score < n3
      #   n4 <= score < n4
      #   score >= n4
      # usage:
      # >> color_ranges [0.3, 0.45, 0.6, 0.88]

    set_formula "a function of a, b, @ka, @kb"
      # setter for the scoring formula
      # expects a string; variables and coefficient: a, @ka; b, @kb
      # usage:
      # >> setformula "1 - (a * @ka) - (b * @kb)"  

    formula?
      # print the scoring formula
      # usage:
      # >> formula?


  Utilities
  ---------
    clear
      # clear the console screen
      # usage:
      # >> clear

    test
      # perform a test scan 
      # usage:
      # >> test
    
EOF
  print manual
  return :done
end

clear