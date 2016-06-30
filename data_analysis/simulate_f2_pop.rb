#encoding: utf-8
require 'bio'
require 'bio-samtools'
require 'pickup'
require 'rinruby'
require 'yaml'
require_relative 'methods_simulate_f2'

if ARGV.empty?
   puts "Please provide directory path of configs.yaml file as argument"
else
   indir = File.expand_path ARGV[0] # location of config file about recombination frequency and number fo chromosomes
end

pars = YAML.load_file("#{indir}/configs.yml")
in_vcf = File.expand_path pars['in_vcf']
xover_file = File.expand_path pars['xovers']
progeny = pars['progeny']
chrs = pars['chrs']
recomb_rate = 0.3

markers = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) } # a hash of variants from vcf file
File.open(in_vcf, 'r').each do |line|
   next if line =~ /^#/
   v = Bio::DB::Vcf.new(line)
   info = v.info
   if info["HET"] == "1"
      markers[v.chrom][v.pos][:ref] = v.ref
      markers[v.chrom][v.pos][:alt] = v.alt
      markers[v.chrom][v.pos][:type] = 'het'
   elsif info["HOM"] == "1"
      markers[v.chrom][v.pos][:ref] = v.ref
      markers[v.chrom][v.pos][:alt] = v.alt
      markers[v.chrom][v.pos][:type] = 'hom'
   end
end

xovers = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) } # a hash of cross over position and prop
File.open(xover_file, 'r').each do |line|
  info = line.split(/\t/)
  next if info[1] !~ /^\d/
  xovers[info[0]][info[1].to_f.ceil] = info[2].to_i
end

# get recombination events in progeny and gametes
chrs = recombinant_progeny(chrs, progeny)
xovers = prop_to_counts(xovers)

gametes = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) } # a hash for recombined gamets

chrs.each_key do | chr |
  chrs[chr][:gametes].each do | num_xos |
    warn "#{chr}\t#{xovers[chr]}\n#{num_xos}"
    recom_pos = recombination_positions(xovers[chr], num_xos.to_i)
    warn "#{recom_pos}"
    gametes = recombined_chromosome(recom_pos, markers[chr])
    warn "#{gametes.class}"
  end
end


