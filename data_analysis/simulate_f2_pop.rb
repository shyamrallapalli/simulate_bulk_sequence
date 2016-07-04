#encoding: utf-8
require 'bio'
require 'bio-samtools'
require 'pickup'
require 'yaml'
require_relative 'methods_simulate_f2'

if ARGV.empty?
   puts "Please provide directory path of configs.yaml file as argument"
else
   indir = File.expand_path ARGV[0] # location of config file about recombination frequency and number fo chromosomes
   Dir.chdir(indir)
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
   markers[v.chrom][v.pos][:ref] = v.ref
   markers[v.chrom][v.pos][:alt] = v.alt
   # info = v.info
   # if info["HET"] == "1"
   #    markers[v.chrom][v.pos][:type] = 'het'
   # elsif info["HOM"] == "1"
   #    markers[v.chrom][v.pos][:type] = 'hom'
   # end
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

# At least in arabidopsis there are more recombinations in male chromosomes than female during meiosis
# http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1002354
sex_recomb_hash = {:male => 6, :female => 4}

gametes = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) } # a hash for recombined gamets

counter = 0 # a counter to index recombined chromosomes at different recombinations per chromosome
chrs.each_key do | chr |
  chrs[chr][:gametes].each do | num_xos |
    sex_recomb_array = []
    sex_recomb_array << Pickup.new(sex_recomb_hash).pick(num_xos.to_i)
    sex_recomb_array.flatten!
    sex_recomb_array.uniq.each do | type |
      count = sex_recomb_array.count(type)
      recom_pos = recombination_positions(xovers[chr], count)
      gametes[chr][count][type][counter] = recombined_chromosome(recom_pos, markers[chr])
    end
    counter += 1
  end
end


