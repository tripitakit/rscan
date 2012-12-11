# signs.rb
# Signs class for rscan program
# author email: patrick.demarta@gmail.com

require 'rubygems'
require 'bio'
require 'rainbow'

class Signs 
  # A multiple-sequence alignment color-masker for signature-sequence search.
  # Scans and scores nucleotide position in multiple sequence alignment input file in fasta format 
  # according to custom grouping, scoring parameters values and formula,  with five color-masking scores range.
  # Groups are defined as an array of arrays of aligned sequences' indexes (admitting Range entries),
  # if omitted at initialize, groups are defined as an array of arrays of each single sequence indexes
  # Scoring parameters are: 
  # scoring_formula :default = 1-(ka/5)*(1-a) - (kb/10)*b
  # a: frequency of nt_base in the ingroup (consenus)
  # b: frequency of nt_base in the outgroup  (aspecificity)
  # ka : consensus coefficient : 20 :strict, 10 :high (majority rule high), 3 :low (majority rule low)
  # kb : aspecificity tollerance coefficient : 0 :allow , 10 :penalty (allowed with penality), 20 :forbid (not allowed)  
  
  attr_reader :input_filename, :alignment, :num_of_seqs, :seq_size, :scores
  attr_accessor :groups, :ka, :kb, :scoring_formula, :score_color_ranges, :paging_window_lenght, :label_size

  def initialize(infile, *groups)
  # required: fasta format multiple sequence alignment filename  
  # optional: groups definition as comma separated list of arrays and/or ranges
    @input_filename = infile    
    @alignment = create_alignment
    @num_of_seqs = @alignment.size
    @seq_size = @alignment[0].seq.size
    @groups = sanitize(groups)
    @ka = 20; @kb = 20
    @scoring_formula = "1 - (@ka*0.5)*(1-a) - (@kb*0.1)*b"
    @score_color_ranges = [0.5, 0.7, 0.8, 0.9]
    @paging_window_lenght = 80
    @label_size = 20
  end


  def create_alignment
  # import a fasta multiple sequences file into an array of Bio::Sequence objects
    alignment = Array.new
    file = Bio::FlatFile.auto(@input_filename)
    file.each do |entry|
      alignment << entry
    end
    return alignment
  end
  
  def sanitize(groups)
    # generate an array of sequence-indexes' arrays,
    # puts each sequence index in a single group-array if groups is empty 
    sanitized_groups = []
    if groups.size  == 0
      @alignment.size.times { |seq_index| sanitized_groups << [seq_index]}
    else
      groups.each { |group|  sanitized_groups <<  group.to_a }
    end
    return sanitized_groups
  end
  
  def scan
  # Iterate along nt positions, calculate nt frequencies for ingroup and outgroups, calculate nt score
    reset_scores!
    sequence_length_range = (0..@seq_size-1)      
    sequence_length_range.each do |position| 
      @groups.size.times do         
        ingroup = @groups.shift.to_a
        outgroup = []
        @groups.each { |group| outgroup << group.to_a }
        outgroup.flatten!
        @groups.push(ingroup)
        ingroup_frequencies = calculate_nt_freq(ingroup, position)
        outgroup_frequencies = calculate_nt_freq(outgroup, position)
        calculate_score(ingroup, position, ingroup_frequencies, outgroup_frequencies)
      end
    end
    return @scores
  end
  
  def reset_scores!
  # create a new scores array of num_of_seqs empty arrays 
    @scores = Array.new
    @num_of_seqs.times { @scores << []}
  end
  
  def calculate_nt_freq(group, position)
  # calculate (a,c,t,g,-) frequencies at current position in group
    num_of_seqs_in_group = group.size
    a = c = t = g = gap = 0
    group.each do |seq_id|
      case base(seq_id, position)
        when "A" then a +=1
        when "C" then c +=1
        when "T" then t +=1
        when "G" then g +=1
        when "-" then gap +=1
      end
    end
    nt_frequencies = {
      "fA" => a.to_f / num_of_seqs_in_group,
      "fT" => t.to_f / num_of_seqs_in_group,
      "fC" => c.to_f / num_of_seqs_in_group,
      "fG" => g.to_f / num_of_seqs_in_group,
      'fgap' => gap.to_f / num_of_seqs_in_group
    }
  end

  def calculate_score(ingroup, position, ingroup_frequencies, outgroup_frequencies)
  # set the current base frequences values on a (consensus) and b (aspecificity) 
  # evaluate the scoring formula and stores the score 
    consensus = specificity = score = 0
    ingroup.each do |seq_id| 
       case base(seq_id, position)
          when "A"
            a = ingroup_frequencies["fA"]
            b = outgroup_frequencies["fA"]
          when "C"
            a = ingroup_frequencies["fC"]
            b = outgroup_frequencies["fC"]
          when "T"
            a = ingroup_frequencies["fT"]
            b = outgroup_frequencies["fT"]
          when "G"
            a = ingroup_frequencies["fG"]
            b = outgroup_frequencies["fG"]
          when "-"
            a = ingroup_frequencies["fgap"]
            b = outgroup_frequencies["fgap"]
        end
      score = eval @scoring_formula
      @scores[seq_id] << score
    end
  end
  
  def color_mask_score(base, nt_score)
  # color maskin according to the base score
    if (nt_score < @score_color_ranges[0])
       base.color(:blue)
    elsif (nt_score >= @score_color_ranges[0] && nt_score < @score_color_ranges[1])
       base.color(:white).background(:blue)
    elsif (nt_score >= @score_color_ranges[1] && nt_score < @score_color_ranges[2])
       base.color(:white).background(:green)
    elsif (nt_score >= @score_color_ranges[2] && nt_score < @score_color_ranges[3])
       base.color(:white).background(:yellow)
    elsif (nt_score >= @score_color_ranges[3] )
       base.color(:white).background(:red)
    end
  end
    
  def color_ruler
    color_ruler = "| < #{@score_color_ranges[0]} ".color(:blue)
    color_ruler += "| .. #{@score_color_ranges[1]} ".color(:white).background(:blue)
    color_ruler += "| .. #{@score_color_ranges[2]} ".color(:white).background(:green)
    color_ruler += "| .. #{@score_color_ranges[3]} ".color(:white).background(:yellow)
    color_ruler += "| > #{@score_color_ranges[3]} |".color(:white).background(:red)
  end
    
  def print_scores
  # prints a paged alignement, with groups separation and color masking
    print "#{color_ruler}\n"
    
    paged_aln = paging
    paged_aln.each_with_index do |alignment_window, page_number|
      print "Page ##{page_number+1}\n"
      @groups.each do |group|
        print page_ruler(page_number)
        group.each do |seq_id|
          paged_seq = alignment_window[seq_id].seq
          print label(seq_id)
          paged_seq.size.times do |position|
            position += @paging_window_lenght * page_number
            base = base(seq_id, position)
            nt_score = @scores[seq_id][position]
            print color_mask_score(base, nt_score)
          end
          puts "\n"
        end
      end  
      puts "\n"
    end
    return true
  end
  
  def paging
  # set the paging windows for the alignment
    num_of_windows = @seq_size / @paging_window_lenght + 1 
    paged_aln =  Array.new
    num_of_windows.times do |window|
      alignment_window = Array.new
      @alignment.each_with_index do |sequence, seq_id|
       window_start = @paging_window_lenght * window
       window_stop = (@paging_window_lenght * window) + @paging_window_lenght - 1 
       window_stop = @seq_size if window_stop > @seq_size
       alignment_window << sequence.seq[window_start..window_stop]
      end
      paged_aln << alignment_window
    end
    return paged_aln
  end
  
  
  def page_ruler(page_number)
  # create the sequence ruler for the current window-page
    ruler = " " * @label_size
    pos = 1
    @paging_window_lenght.times do
      if pos/10.0 == pos/10 
        ruler += "|" 
      else
        ruler += "-"
      end
      pos += 1
    end
    ruler += (@paging_window_lenght + page_number*@paging_window_lenght).to_s + "\n"
    ruler.color(:magenta)
  end
      
  def label(seq_id)
  # format seq labels to @label_size chars lenght
    definition = seq_id.to_s + ". " + @alignment[seq_id].definition.to_s
    label = if (definition.size <= @label_size)
      definition + " " * (@label_size-definition.size)
    else
      definition[0..@label_size] 
    end
    label.foreground(:blue)
  end

  def base(seq_id, position)
    # retun the uppercase base given seq_id and position
    @alignment[seq_id].seq[position,1].upcase
  end

  def scores2csv(outfile)
  # write the scores in a csv outfile 
    CSV.open(outfile, "wb") do |csv|
      @scores.each_with_index do |seq_score, seq_index|
        definition = @alignment[seq_index].definition
        csv << [definition].concat(seq_score)
      end
    end
  end

  def groups
  # prints the groups with sequences' labels
    @groups.each_with_index do |group, grp_index|
      puts "Group #{grp_index}: #{group.inspect}".color(:blue)
      group.to_a.each do |seq_index|
        puts "\t" + seq_index.to_s + ". " + @alignment[seq_index].definition
      end
    end
  end

  def labels
  # prints only the sequences' labels 
    @alignment.each_with_index do |sequence, index|
      puts index.to_s + ". " + sequence.definition
    end
    return true
  end
  
end