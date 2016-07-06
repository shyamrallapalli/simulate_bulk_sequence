#encoding: utf-8
require 'bio'
require 'bio-samtools'


if ARGV.empty?
   puts "Please provide a fasta file, a vcf file and an short phrase to include in filename as arguments in that order"
   puts "usage: ruby fasta_var_udpate.rb fasta vcf sample"
   exit
else
   in_fasta = File.expand_path ARGV[0]
   in_vcf = File.expand_path ARGV[1]
   sample = ARGV[2].chomp
end

# a hash of sequences from fasta file
sequences = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }
Bio::FastaFormat.open(in_fasta).each do |fas|
  fas.definition += ' ' + sample
  sequences[fas.entry_id][:def] = fas.definition
  sequences[fas.entry_id][:seq] = fas.seq
end

# a hash of variants from vcf file
variants = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }
File.open(in_vcf, 'r').each do |line|
   next if line =~ /^#/
   v = Bio::DB::Vcf.new(line)
   variants[v.chrom][v.pos][:ref] = v.ref
   variants[v.chrom][v.pos][:alt] = v.alt
end

indels = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }
variants.each_key do | chr |
  variants[chr].each_key do | pos |
    ref = variants[chr][pos][:ref]
    alt = variants[chr][pos][:alt]
    if ref.length == alt.length
      # string index starts at '0', while positions start at '1'
      sequences[chr][:seq][pos-1] = alt
    else # indels
      indels[chr][pos] = variants[chr][pos]
    end
  end
end

indels.each_key do | chr |
  # decreasing order of positions
  sorted_pos = indels[chr].keys.sort { |a, b| b <=> a }
  sorted_pos.each do | pos |
    len = indels[chr][pos].length
    if length == 1
      sequences[chr][:seq][pos-1] = alt
    else
      sequences[chr][:seq][pos-1..len-1] = alt
    end
  end
end

out_fasta = File.open(in_fasta + '_' + sample + '.fas', 'w')
sequences.each_key do | chr |
  seq = sequences[chr][:seq]
  id = sequences[chr][:def]
  out_fasta.puts seq.to_fasta(id, 80)
end
out_fasta.close

