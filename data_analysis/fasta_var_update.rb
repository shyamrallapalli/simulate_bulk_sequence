#encoding: utf-8
require 'bio'
require 'bio-samtools'


def write_variant_to_chr(variants, fas_entry, outfile)
  chr = fas_entry.entry_id
  indels = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }
  variants[chr].each_key do | pos |
    ref = variants[chr][pos][:ref]
    alt = variants[chr][pos][:alt]
    if ref.length == alt.length
      # string index starts at '0', while positions start at '1'
      fas_entry.seq[pos-1] = alt
    else # indels
      indels[pos] = variants[chr][pos]
    end
  end

  # decreasing order of positions
  sorted_pos = indels.keys.sort { |a, b| b <=> a }
  sorted_pos.each do | pos |
    len = indels[pos].length
    if length == 1
      fas_entry.seq[pos-1] = alt
    else
      stop = spos - 1 - len - 1
      fas_entry.seq[pos-1..stop] = alt
    end
  end

  outfile.puts fas_entry.seq.to_fasta(fas_entry.definition, 79)
end

if ARGV.empty?
   puts "Please provide a vcf file, a fasta file and an short phrase to include in new filename as arguments in that order"
   puts "usage: ruby fasta_var_udpate.rb vcf fasta sample"
   exit
else
   in_vcf = File.expand_path ARGV[0]
   in_fasta = File.expand_path ARGV[1]
   sample = ARGV[2].chomp
end

# a hash of variants from vcf file
variants = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }
File.open(in_vcf, 'r').each do |line|
   next if line =~ /^#/
   v = Bio::DB::Vcf.new(line)
   variants[v.chrom][v.pos][:ref] = v.ref
   variants[v.chrom][v.pos][:alt] = v.alt
end

out_fasta = File.open(in_fasta + '_' + sample + '.fas', 'w')
# a hash of sequences from fasta file
sequences = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }
Bio::FastaFormat.open(in_fasta).each do |fas|
  fas.definition += ' ' + sample
  write_variant_to_chr(variants, fas, out_fasta)
end
out_fasta.close
