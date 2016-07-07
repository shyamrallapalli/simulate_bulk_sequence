#encoding: utf-8
require 'bio'
require 'bio-samtools'
require_relative 'update_chr_seq'

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
  fas = update_variant_to_chr(variants, fas)
  out_fasta.puts fas.seq.to_fasta(fas.definition, 79)
end
out_fasta.close
