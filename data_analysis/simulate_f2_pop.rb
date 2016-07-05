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
progeny_num = pars['progeny']
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
chrs = recombinant_progeny(chrs, progeny_num)
xovers = prop_to_counts(xovers)

gametes = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) } # a hash for recombined gamets

counter = 0 # a counter to index recombined chromosomes at different recombinations per chromosome
chrs.each_key do | chr |
  chrs[chr][:gametes].each do | num_xos |
    gender_recomb_hash = recombinant_gender_num(num_xos)
    gender_recomb_hash.each_key do | type |
      count = gender_recomb_hash[type]
      recom_pos = recombination_positions(xovers[chr], count)
      gametes[chr][count][type][counter] = recombined_chromosome(recom_pos, markers[chr])
    end
    counter += 1
  end
end

# warn "#{gametes}"

def get_male_female_index(gametes, chr, num_xos)
  male = gametes[chr][num_xos][:male].keys.sample
  female = gametes[chr][num_xos][:female].keys.sample
  until male != female
    female = gametes[chr][num_xos][:female].keys.sample
  end
  [male,female]
end

progeny = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) } # a hash for recombined gamets
counter = 0
chrs.each_key do | chr |
  chrs[chr][:progeny].each do | num_xos |
    if num_xos == 0
      one, two = randomize_pair
      male, female = get_male_female_index(gametes, chr, num_xos)
      progeny[chr][counter][:male] = gametes[chr][num_xos][:male][male][one]
      gametes[chr][num_xos][:male][male].delete(one)
      progeny[chr][counter][:female] = gametes[chr][num_xos][:female][female][two]
      gametes[chr][num_xos][:female][female].delete(two)
    else
      gender_recomb_array = recombinant_gender_num(num_xos)
      one, two = randomize_pair
      gender_recomb_array.uniq.each do | type |
        count = gender_recomb_array.count(type)
        index = gametes[chr][count][type].keys.sample
        progeny[chr][counter][type] = gametes[chr][count][type][index][one]
        gametes[chr][count][type][index].delete(one)
        one = two
      end
    end
  end
  counter += 1
end

warn "#{progeny}"


